from flask import Flask, redirect, url_for, request
app = Flask(__name__)


################ add the different routes @###################

@app.route('/')
def hello(): #here we can name the function what we want as long as the app.route call decorator is correctly called by flask
 return 'Hello from the main'

@app.route('/hello')
def hello_world():
 return 'Hello World'

@app.route('/hello/<name>')
def hello_name(name):
 return 'Hello '+str(name)

@app.route('/timesTwo/<int:postID>')
def timesTwo(postID):
 return 'Result '+str(postID*2)

@app.route('/timesFour/<int:po>')
def timesFour(po):
 return 'Result '+str(po*4)

#################### GET/POST ##############
@app.route('/success/<name>')
def success(name):
 return 'welcome %s' % name

@app.route('/login',methods = ['POST', 'GET'])
def login():
 if request.method == 'POST':
  user = request.form['nm']
  return redirect(url_for('success',name = user))
 else:
  user = request.args.get('nm')
  return redirect(url_for('success',name = user))

if __name__ == '__main__':
 app.run(debug=True)
