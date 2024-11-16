FROM picker24/root_v6_26_10:alma9 AS nuwro_build
WORKDIR /

ARG NUWRO_VERSION=master
ENV NUWRO_VERSION=${NUWRO_VERSION}

ENV NUWRO_WORK_DIR=/opt/nuwro-src/

WORKDIR /opt/
RUN --mount=type=ssh git clone https://github.com/NuWro/nuwro.git ${NUWRO_WORK_DIR}
WORKDIR /opt/nuwro-src
RUN git checkout ${NUWRO_VERSION}

RUN sed -i "s:-lPythia6:-L${PYTHIA6} -lPythia6:g" Makefile
RUN make -j4

FROM picker24/root_v6_26_10:alma9
WORKDIR /

RUN mkdir -p /opt/nuwro/build/src/dis
RUN mkdir -p /opt/nuwro/build/src/espp
RUN mkdir -p /opt/nuwro/build/src/sf
RUN mkdir -p /opt/nuwro/build/src/rpa
RUN mkdir -p /opt/nuwro/build/src/rew
COPY --from=nuwro_build /opt/nuwro-src/data /opt/nuwro/build/data
COPY --from=nuwro_build /opt/nuwro-src/bin /opt/nuwro/build/bin
COPY --from=nuwro_build /opt/nuwro-src/src/*.h /opt/nuwro/build/src/
COPY --from=nuwro_build /opt/nuwro-src/src/dis/*.h /opt/nuwro/build/src/dis/
COPY --from=nuwro_build /opt/nuwro-src/src/espp/*.h /opt/nuwro/build/src/espp/
COPY --from=nuwro_build /opt/nuwro-src/src/sf/*.h /opt/nuwro/build/src/sf/
COPY --from=nuwro_build /opt/nuwro-src/src/rpa/*.h /opt/nuwro/build/src/rpa/
COPY --from=nuwro_build /opt/nuwro-src/src/rew/*.h /opt/nuwro/build/src/rew/

ENV NUWRO=/opt/nuwro/build
ENV PATH=${NUWRO}/bin:${PATH}
ENV ROOT_INCLUDE_PATH=${NUWRO}/src:${ROOT_INCLUDE_PATH}
