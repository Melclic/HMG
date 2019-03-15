An overly complicated webservice for the libCH package.

TODO:
Ultimetely would like to have everything in a docker (including RD-kit)

docker build -t testch .
docker run --name testch -p 5000:5000 --rm testch 

curl -H "Content-type: application/json" -X POST http://localhost:5000/REST -d '{"tau": "27.0", "C": "40.0", "D": "20.0", "a": "0.3"}'
