#This pulls an already existing docker container (default here is the docker image from the web for python 3.6 with alpine linux that is a smaller version of linux compared to more heavy ones as ubuntu)
FROM python:3.6-alpine

#add a user (maybe send the command to create a group and add that user to the group)
RUN adduser -D webCH

WORKDIR /home/webCH
COPY . /home/webCH

#install some of the required packages
#alternatively can have all the files installed in a requirements.txt file
#				--> should be frozen from a working version of the code
RUN pip install flask
RUN pip install flask-restful
#should have conda 

#This is for testing that the docker does not actually have root access
#CMD ["cat", "/tmp/secrets.txt"]

#at this point all the files are still root. This makes sure that the permission granted to the whole files in the container are webCH user
RUN chown -R webCH:webCH ./
USER webCH

#Make the docker
EXPOSE 5000

#run the python script for the webservice
#NOTE this should be a bash executable file that contains permision
#ENTRYPOINT ["python"]

CMD ["python", "app.py"]
