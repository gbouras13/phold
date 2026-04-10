
FROM quay.io/pawsey/pytorch:2.7.1-rocm6.3.3

## install foldseek only non python dep

RUN apt-get update && apt-get install -y \
    wget \
    tar \
    aria2 \
    && rm -rf /var/lib/apt/lists/*

# install foldseek
RUN wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && \
    tar -xzf foldseek-linux-avx2.tar.gz && \
    rm foldseek-linux-avx2.tar.gz && \
    mv foldseek /opt/foldseek

ENV PATH="/opt/foldseek/bin:${PATH}"

# install phold

RUN pip install --ignore-installed PyYAML
RUN python3 -c "import torch; print(torch.__version__)" || echo "PyTorch is not installed before installing the project. "
ARG PHOLD_REF
RUN pip install git+https://github.com/gbouras13/phold.git@${PHOLD_REF}
RUN pip install pytest

# tests
RUN foldseek --help
RUN phold -h
RUN phold run -h
