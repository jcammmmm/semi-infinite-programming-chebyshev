import markdown as md
from jinja2 import Environment, FileSystemLoader

env = Environment(
  loader=FileSystemLoader('.'),
  autoescape=False
)

with open('content.md', 'r', encoding='UTF-8') as input_file:
  text = input_file.read()
html = md.markdown(text, output_format='html5', tab_length=2)
print(html)

post = env.get_template('layout.html').render(content=html)
open('content.html', 'w+').write(post)
