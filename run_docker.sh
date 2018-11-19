#!/bin/bash
docker run -it -p 8888:8888 \
    -v $(pwd):/home/jovyan \
    -v $(pwd)/data:/home/jovyan/data \
    jupyter/scipy-notebook:e8613d84128b

    