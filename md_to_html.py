import markdown as md
from jinja2 import Environment, FileSystemLoader

# FILENAME = '220808Fileshare'
FILENAME = 'content'
env = Environment(
  loader=FileSystemLoader('.'),
  autoescape=False
)

with open(FILENAME + '.md', 'r', encoding='UTF-8') as input_file:
  text = input_file.read()
html = md.markdown(text, output_format='html5', tab_length=2)
print(html)

post = env.get_template('layout.html').render(content=html)
open(FILENAME + '.html', 'w+').write(post)