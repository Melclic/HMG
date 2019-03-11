

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

	###### from job_conf.xml.sample_advanced ########	

    <destinations default="local">
        <!-- Destinations define details about remote resources and how jobs
             should be executed on those remote resources.
         -->
        <destination id="local" runner="local"/>
        <destination id="multicore_local" runner="local">
          <param id="local_slots">4</param> <!-- Specify GALAXY_SLOTS for local jobs. -->
          <!-- Warning: Local slot count doesn't tie up additional worker threads, to prevent over
               allocating machine define a second local runner with different name and fewer workers
               to run this destination. -->
          <param id="embed_metadata_in_job">True</param>
          <!-- Above parameter will be default (with no option to set
               to False) in an upcoming release of Galaxy, but you can
               try it early - it will slightly speed up local jobs by
               embedding metadata calculation in job script itself.
          -->
          <job_metrics />
          <!-- Above element demonstrates embedded job metrics definition - see
               job_metrics_conf.xml.sample for full documentation on possible nested
               elements. This empty block will simply disable job metrics for the
               corresponding destination. -->
        </destination>
        <destination id="docker_local" runner="local">
          <param id="docker_enabled">true</param>
          <!-- docker_volumes can be used to configure volumes to expose to docker,
               For added isolation append :ro to the path to mount it read only.
               Galaxy will attempt to infer a reasonable set of defaults which
               volumes should be exposed how based on Galaxy's settings and the
               destination - but be sure to add any library paths or data incides
               that may be needed read-only.
          -->
          <!--
          <param id="docker_volumes">$defaults,/mnt/galaxyData/libraries:ro,/mnt/galaxyData/indices:ro</param>
          -->
          <!-- For a stock Galaxy instance and traditional job runner $defaults will
               expand out as:

               $galaxy_root:ro,$tool_directory:ro,$job_directory:ro,$working_directory:rw,$default_file_path:rw

               This assumes most of what is needed is available under Galaxy's root directory,
               the tool directory, and the Galaxy's file_path (if using object store creatively
               you will definitely need to expand defaults).

               This configuration allows any docker instance to write to any Galaxy
               file - for greater isolation set outputs_to_working_directory in
               galaxy.ini. This will cause $defaults to allow writing to much
               less. It will then expand as follows:

               $galaxy_root:ro,$tool_directory:ro,$job_directory:ro,$working_directory:rw,$default_file_path:ro

               If using the Pulsar, defaults will be even further restricted because the
               Pulsar will (by default) stage all needed inputs into the job's job_directory
               (so there is not need to allow the docker container to read all the
               files - let alone write over them). Defaults in this case becomes:

               $job_directory:ro,$tool_directory:ro,$job_directory/outputs:rw,$working_directory:rw

               Python string.Template is used to expand volumes and values $defaults,
               $galaxy_root, $default_file_path, $tool_directory, $working_directory,
               are available to all jobs and $job_directory is also available for
               Pulsar jobs.
          -->
          <!-- One can run docker using volumes-from tag by setting the following
               parameter. For more information on volumes-from check out the following
               docker tutorial. https://docs.docker.com/userguide/dockervolumes/
          -->
          <!-- <param id="docker_volumes_from">parent_container_name</param> -->
          <!-- Control memory allocatable by docker container with following option:
          -->
          <!-- <param id="docker_memory">24G</param> -->
          <!-- By default Docker will need to runnable by Galaxy using
               password-less sudo - this can be configured by adding the
               following line to the sudoers file of all compute nodes
               with docker enabled:

               galaxy  ALL = (root) NOPASSWD: SETENV: /usr/bin/docker

               The follow option is set to false to disable sudo (docker
               must likewise be configured to allow this).
          -->
          <!-- <param id="docker_sudo">false</param> -->
          <!-- Following option can be used to tweak sudo command used by
               default. -->
          <!-- <param id="docker_sudo_cmd">/usr/bin/sudo -extra_param</param> -->
          <!-- By default, docker container will not have any networking
               enabled. host networking can be bridged by uncommenting next option
               http://docs.docker.io/reference/run/#network-settings
          -->
          <!-- <param id="docker_net">bridge</param> -->
          <!-- By default, a container will live on past its run. By
               adding the '\-\-rm' flag to the command line, the container
               will be removed automatically after the program is complete.
          -->
          <!-- <param id="docker_auto_rm">true</param> -->
          <!-- Override which user to launch Docker container as - defaults to 
               Galaxy's user id. For remote job execution (e.g. Pulsar) set to
               remote job user. Leave empty to not use the -u argument with
               Docker. -->
          <!-- <param id="docker_set_user">$UID</param> -->
          <!-- Pass extra arguments to the docker run command not covered by the
               above options. -->
          <!-- <param id="docker_run_extra_arguments"></param> -->
          <!-- Following command can be used to tweak docker command. -->
          <!-- <param id="docker_cmd">/usr/local/custom_docker/docker</param> -->
          <!-- Following can be used to connect to docke server in different
               ways (translated as -H argument to docker client). -->
          <!-- <param id="docker_host">unix:///var/run/docker.sock</param> -->
          <!-- <param id="docker_host">:5555</param> -->
          <!-- <param id="docker_host">:5555</param> -->
          <!-- <param id="docker_host">tcp://127.0.0.1:4243</param> -->

          <!-- If deployer wants to use docker for isolation, but does not
               trust tool's specified container - a destination wide override
               can be set. This will cause all jobs on this destination to use
               that docker image. -->
          <!-- <param id="docker_container_id_override">busybox:ubuntu-14.04</param> -->

          <!-- Likewise, if deployer wants to use docker for isolation and
               does trust tool's specified container - but also wants tool's not
               configured to run in a container the following option can provide
               a fallback. -->
          <!-- <param id="docker_default_container_id">busybox:ubuntu-14.04</param> -->

          <!-- If the destination should be secured to only allow containerized jobs
               the following parameter may be set for the job destination. Not all,
               or even most, tools available in Galaxy core or in the Tool Shed
               support Docker yet so this option may require a lot of extra work for
               the deployer. -->
          <!-- <param id="require_container">true</param> -->
        </destination>
