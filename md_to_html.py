import markdown as md

with open('doc.md', 'r', encoding='UTF-8') as input_file:
  text = input_file.read()
html = md.markdown(text, output_format='html5', tab_length=2)
print(html)

with open('doc.html', 'w', encoding='UTF-8') as output_file:
  output_file.write(html)