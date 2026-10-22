# Reproducible container. Pin the Bioconductor release you validated with.
FROM bioconductor/bioconductor_docker:RELEASE_3_21

WORKDIR /app
COPY install.R /app/install.R
RUN Rscript install.R

COPY . /app
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"]
