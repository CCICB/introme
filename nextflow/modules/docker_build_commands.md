
```
docker build -f nextflow/modules/pangolin/Dockerfile.cpu  \
    --build-arg IMAGE_VERSION=2.0.0 \
    -t introme-pangolin:2.0.0-runtime-cpu     .

docker build -f nextflow/modules/pangolin/Dockerfile.gpu \
    --build-arg IMAGE_VERSION=2.0.0 \
    -t introme-pangolin:2.0.0-runtime-gpu     .

docker build -f nextflow/modules/pangolin/Dockerfile.annotation \
    --build-arg BASE_IMAGE=introme-pangolin:2.0.0-runtime-cpu \
    --build-arg GENCODE_VERSION=50 \
    --build-arg GENCODE_GTF=</path/to/gtf>
    -t introme-pangolin:2.0.0-cpu-gencode50     .

```
