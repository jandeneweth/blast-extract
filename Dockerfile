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
# DEBUG: mount the scripts directory for latest versions!
#COPY ./scripts ./scripts


ENTRYPOINT []
CMD []
