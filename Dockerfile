FROM jupyter/scipy-notebook:e8613d84128b

USER root

RUN pip install cmapPy==3.3.3 \
	umap==0.1.1


