FROM ipython/scipyserver

RUN python -m pip install Cython==0.23 
RUN python -m pip install cmapPy==2.2.2 \
	umap==0.1.1 \
	seaborn==0.9.0 \
	joblib==0.13.0

