

## GALAXY TOOL

1) build the image using the DockerFile

  docker build -t ibisba/ch .

2) Add tool to tool_conf.xml

  <section id="retro" name="Retro Nodes">
    <tool file="/home/mdulac/workspace/HMG/libCH/ch_docker_galaxy.xml" />
  </section>

## GALAXY DOCKER

1) Added to tool xml:

  <requirements>
    <container type="docker">ibibsa/ch:latest</container>
  </requirements>

2) Need to create/modify config/job_conf.xml:
	
  <destinations default="docker_local">
    <destination id="local" runner="local" />
    <destination id="docker_local" runner="local">
      <param id="docker_enabled">true</param>
      <param id="docker_sudo">false</param> #<--- this is assuming that docker can be run without a password
    </destination>
  </destinations>


3) to the python script add: #!/usr/bin/env python3 

## PULSAR






#Draft

Here we have:
	- Tool running inside a docker
	- Flask server that accepts JSON REST requests
	- A job handling (celery?, REDIS queue?) TODO
Here we have the Cooper-Helmstetter cell cycle model

gcc -shared -o cooperHelmstetter.so -fPIC ori.c -lm

docker build -t ibisba/ch .
docker run --name testch -p 5000:5000 --rm testch

curl -H "Content-type: application/json" -X POST http://localhost:5000/REST -d '{"tau": "27.0", "C": "40.0", "D": "20.0", "a": "0.3"}'



2) Inside galaxy's config:

  cp job_conf.xml.sample_basic job_conf.xml

  <destinations default="docker_local">
       <destination id="local" runner="local"/>
       <destination id="docker_local" runner="local">
         <param id="docker_enabled">true</param>
       </destination>
    </destinations>
