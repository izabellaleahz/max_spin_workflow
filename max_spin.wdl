version 1.0
workflow max_spin {
    input {
        String output_directory
        File anndata_file
        String sample_key = 'sample'

        #general parameters
        Int cpu = 24
        String memory = "128G"
        Int extra_disk_space = 32
        String docker = "mparikhbroad/hotspot:latest"
        Int preemptible = 2
    }
    String output_directory_stripped = sub(output_directory, "/+$", "")
    call run_max_spin {
        input:
            output_dir = output_directory_stripped,
            anndata_file = anndata_file,
            sample_key = sample_key,
            
            cpu=cpu,
            memory=memory,
            extra_disk_space = extra_disk_space,
            docker=docker,
            preemptible=preemptible
    }
    output {
        File maxspin_object = run_max_spin.maxspin_object
    }
}
task run_max_spin {
    input {
        String output_dir
        File anndata_file
        String sample_key

        String memory
        Int extra_disk_space
        Int cpu
        String docker
        Int preemptible
    }
    command <<<
        set -e
        mkdir -p outputs
        python <<CODE

        #imports 
        from maxspin import spatial_information, pairwise_spatial_information
        import numpy as np
        import scanpy as sc
        import squidpy as sq
        import scvi

        #load in anndata object 
        adata_full = sc.read_h5ad('~{anndata_file}')

        samples = list(set(adata.obs['~{sample_key}']))

        for sample in sample: 
            adata_list.append(adata_full[adata_full.obs['~{sample_key}'] == sample])

        # look at autocorrelation per sample 
        
        scvi.model.SCVI.setup_anndata(adata_list)
        model = scvi.model.SCVI(adata_list, n_latent=20)

        # Sample log-expression values from the posterior.
        posterior_samples = np.log(model.get_normalized_expression(return_numpy=True, return_mean=False, n_samples=20, library_size="latent"))
        adata_scvi = adata_list.copy()
        adata_scvi.X = np.mean(posterior_samples, axis=0)
        adata_scvi.layers["std"] = np.std(posterior_samples, axis=0)

        spatial_information(adata_scvi, prior="gaussian")
        #export an anndata object per sample 
        adata_scvi.write_h5ad('output/autocorrelation_output.h5ad')
        
        #pairwise correlation 
        Nnmf = 20
        nmf = NMF(n_components=Nnmf, init='nndsvd', random_state=0)
        W = nmf.fit_transform(np.exp(adata_scvi.X))
        adata_nmf = sc.AnnData(
            X=W,
            layers={"log1p": np.log1p(W)},
            obsm=adata.obsm,
            obs=adata.obs,
            obsp=adata.obsp,
            uns=adata.uns
            var=adata.var)
        #computer spatial information for each meta gene 
        spatial_information(adata_nmf, layer="log1p", prior=None)
        #pairwise spatial information 
        pairwise_spatial_information(adata_nmf, layer="log1p", prior=None)
        #export an anndata object per sample 
        adata_scvi.write_h5ad('output/pairwise_output.h5ad')

        CODE
        gsutil -m rsync -r outputs ~{output_dir}
    >>>
    output {
        File maxspin_object = 'outputs/'
    }
    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + (ceil(size(anndata_file, "GB")*4) + extra_disk_space) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }
}