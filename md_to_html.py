import markdown as md
import datetime as time
from jinja2 import Environment, FileSystemLoader

# FILENAME = '220808Fileshare'
FILENAME = 'chebyshev'
env = Environment(
  loader=FileSystemLoader('.'),
  autoescape=False
)

with open(FILENAME + '.md', 'r', encoding='UTF-8') as input_file:
  text = input_file.read()
'''
attr_list: for inline attribute definitions
toc      : table of contents
tables   : generates html tables
'''
html = md.markdown(text, extensions=['attr_list', 'toc', 'tables'], output_format='html5', tab_length=2, )
print(html)

post = env.get_template('layout.html').render(content=html)
open(FILENAME + '.html', 'w+', encoding='UTF-8').write(post)
open(FILENAME + '.html', 'a').write('last update: ' + str(time.datetime.now()))