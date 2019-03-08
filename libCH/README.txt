Here we have:
	- Tool running inside a docker
	- Flask server that accepts JSON REST requests
	- A job handling (celery?, REDIS queue?) TODO
Here we have the Cooper-Helmstetter cell cycle model

gcc -shared -o cooperHelmstetter.so -fPIC ori.c -lm

docker build -t testch .
docker run --name testch -p 5000:5000 --rm testch

curl -H "Content-type: application/json" -X POST http://localhost:5000/REST -d '{"tau": "27.0", "C": "40.0", "D": "20.0", "a": "0.3"}'
