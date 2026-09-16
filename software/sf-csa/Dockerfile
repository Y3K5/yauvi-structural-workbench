# A runnable SF-CSA image: the package plus its two external binaries, pinned.
#
# Deliberately campaign-free. It builds the component and nothing else — no
# campaign root, no project mounts, no manifests. Those are supplied at run time,
# which is what lets one image serve any campaign.
#
# The version pins matter: `required_foldseek_version` in a database manifest is
# checked at run time and a mismatch raises, because an aligner version change is
# a change to the results.
FROM mambaorg/micromamba:2.0.5

USER root
WORKDIR /opt/sf-csa

COPY pyproject.toml README.md ./
COPY src ./src

RUN micromamba install -y -n base -c conda-forge -c bioconda \
      python=3.13 foldseek=10.941cd33 diamond=2.2.5 && \
    micromamba clean --all --yes && \
    /opt/conda/bin/python -m pip install --no-cache-dir .

ENV PATH=/opt/conda/bin:$PATH
USER mambauser

ENTRYPOINT ["sf-csa"]
CMD ["--help"]
