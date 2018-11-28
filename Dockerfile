FROM jupyter/scipy-notebook:e8613d84128b

RUN python -m pip install umap-learn==0.3.6 \
	joblib==0.13.0 \
	plotly==3.4.2 \
	python-igraph==0.7.1

