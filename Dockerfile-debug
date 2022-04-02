FROM python:3.10.2-bullseye

RUN apt-get update

# Set up working directory
WORKDIR /usr/src/app
RUN mkdir ./tools

# Get NCBI BLAST
RUN apt-get -y install ncbi-blast+

# Install python libraries
COPY requirements.txt ./
RUN python3 -m pip install --no-cache-dir -r requirements.txt
RUN rm requirements.txt

# Add pipeline scripts
# DEBUG: mount a site-packages directory for latest versions without re-building image!
ENV PYTHONPATH "${PYTHONPATH}:/usr/src/app/site-packages"


ENTRYPOINT []
CMD []
