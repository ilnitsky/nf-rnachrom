FROM continuumio/miniconda3:latest

WORKDIR /app

# Copy the environment file
COPY envs/secondary_processing.yml .

# Create the conda environment
RUN conda env create -f secondary_processing.yml

# Make RUN commands use the new environment
SHELL ["conda", "run", "-n", "secondary_processing", "/bin/bash", "-c"]

# Set entrypoint to activate the conda environment
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "secondary_processing"]

# Default command - can be overridden
CMD ["/bin/bash"] 