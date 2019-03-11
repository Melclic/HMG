#This pulls an already existing docker container (default here is the docker image from the web for python 3.6 with alpine linux that is a smaller version of linux compared to more heavy ones as ubuntu)
FROM python:3.6-alpine

#add a user (maybe send the command to create a group and add that user to the group)
#RUN adduser -D restCH


#install some of the required packages
#alternatively can have all the files installed in a requirements.txt file
#				--> should be frozen from a working version of the code
RUN apk update
RUN apk add build-base
#RUN pip install flask
#RUN pip install flask-restful
#should have conda 

#This is for testing that the docker does not actually have root access
#CMD ["cat", "/tmp/secrets.txt"]

#at this point all the files are still root. This makes sure that the permission granted to the whole files in the container are webCH user
#RUN chown -R restCH:restCH ./
#USER restCH

WORKDIR /home/ch
COPY . /home/ch
RUN gcc -shared -o cooperHelmstetter.so -fPIC ori.c -lm

### This simplifies the call to the docker as the scripts are stored in the bin
RUN cp cooperHelmstetter.so /usr/bin/
RUN cp cooperHelmstetter_local.py /usr/bin/
RUN chmod a+x /usr/bin/cooperHelmstetter_local.py
RUN chmod a+x /usr/bin/cooperHelmstetter.so

#Make the docker
#EXPOSE 5000

#run the python script for the webservice
#NOTE this should be a bash executable file that contains permision
#ENTRYPOINT ["python"]

#RUN chown -R $(whoami) /usr/local/bin

# add smalt tools to /usr/bin so they're on the PATH
#RUN ln -s cooperHelmstetter_local.py /usr/bin
#RUN ln -s /tool/smalt_wrapper.py /usr/bin

# make smalt_wrapper executable
#RUN chmod a+x cooperHelmstetter_local.py

#CMD ["python", "cooperHelmstetter_local.py"]
#CMD ["python", "cooperHelmstetter_flask.py"]
#CMD ["/bin/bash"]

#ENTRYPOINT ["sh"]
