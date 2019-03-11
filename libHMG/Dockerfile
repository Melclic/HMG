#This pulls an already existing docker container (default here is the docker image from the web for python 3.6 with alpine linux that is a smaller version of linux compared to more heavy ones as ubuntu)
FROM python:3.6-alpine

#add a user (maybe send the command to create a group and add that user to the group)
#RUN adduser -D flaskHMG


RUN apk update
RUN apk add build-base
RUN apk add gsl-dev
#RUN pip install flask
#RUN pip install flask-restful

#This is for testing that the docker does not actually have root access
#CMD ["cat", "/tmp/secrets.txt"]

#at this point all the files are still root. This makes sure that the permission granted to the whole files in the container are flaskHMG user
#RUN chown -R flaskHMG:flaskHMG ./
#USER flaskHMG

WORKDIR /home/hmg
COPY . /home/hmg
RUN gcc -shared -o abm.so -fPIC pyConnect.c model.c cellPopulation.c cell.c utility.c inputModel.c -lm -lgsl -lgslcblas

### This simplifies the call to the docker as the scripts are stored in the bin
RUN cp abm.so /usr/bin/
RUN cp abm_local.py /usr/bin/
RUN chmod a+x /usr/bin/abm_local.py
RUN chmod a+x /usr/bin/abm.so

#Make the docker
#EXPOSE 5000

#run the python script for the webservice
#NOTE this should be a bash executable file that contains permision
#ENTRYPOINT ["python"]

#CMD ["python", "abm.py"]
#ENTRYPOINT ["bash"]
