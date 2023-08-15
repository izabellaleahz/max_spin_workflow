version 1.0
workflow maxspin {
    input {
        String output_directory
        File anndata_file
        String sample_key = 'sample'

        #general parameters
        Int cpu = 24
        String memory = "128G"
        Int extra_disk_space = 32
        String docker = "us.gcr.io/landerlab-atacseq-200218/maxspin:latest"
        Int preemptible = 2
    }
    String output_directory_stripped = sub(output_directory, "/+$", "")
    call run_maxspin {
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
        File maxspin_object = run_maxspin.maxspin_object
    }
}
task run_maxspin {
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
        from maxspin import spatial_information, pairwise_spatial_information
        import numpy as np
        import scanpy as sc
        import squidpy as sq
        import scvi
        adata_full = sc.read_h5ad('~{anndata_file}')

        scvi.model.SCVI.setup_anndata(adata_full)
        model = scvi.model.SCVI(adata_full, n_latent=20)
        posterior_samples = np.log(model.get_normalized_expression(return_numpy=True, return_mean=False, n_samples=20, library_size="latent"))
        adata_scvi = adata_list.copy()
        adata_scvi.X = np.mean(posterior_samples, axis=0)
        adata_scvi.layers["std"] = np.std(posterior_samples, axis=0)
        spatial_information(adata_scvi, prior="gaussian")
        adata_scvi.write_h5ad('~{output_dir}' + 'autocorrelation_output.h5ad')
        Nnmf = 20
        nmf = NMF(n_components=Nnmf, init='nndsvd', random_state=0)
        W = nmf.fit_transform(np.exp(adata_scvi.X))
        adata_nmf = sc.AnnData(
            X=W,
            layers={"log1p": np.log1p(W)},
            obsm=adata.obsm,
            obs=adata.obs,
            obsp=adata.obsp,
            uns=adata.uns,
            var=adata.var)
        spatial_information(adata_nmf, layer="log1p", prior=None)
        pairwise_spatial_information(adata_nmf, layer="log1p", prior=None)
        adata_scvi.write_h5ad('~{output_dir}' + 'pairwise_output.h5ad')
        CODE
        gsutil -m rsync -r outputs ~{output_dir}
    >>>
    output {
        File maxspin_object = '~{output_dir}' + 'pairwise_output.h5ad'
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