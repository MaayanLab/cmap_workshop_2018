FROM ipython/scipyserver

# RUN python -m pip install Cython==0.23 
RUN python -m pip install umap-learn \
	seaborn==0.9.0 \
	joblib==0.13.0 \
	plotly==3.4.2

