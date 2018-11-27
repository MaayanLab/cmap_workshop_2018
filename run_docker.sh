#!/bin/bash
docker run -it -p 8888:8888 \
	-e "PASSWORD=password" -e "USE_HTTP=1" \
    -v $(pwd):/notebooks \
    -v $(pwd)/data:/notebooks/data \
    maayanlab/ipython2-cmap

    