from html.parser import HTMLParser
from pathlib import Path
import re,json,hashlib,datetime

class Node:
 def __init__(self,tag='',attrs=()):self.tag=tag;self.attrs=dict(attrs);self.children=[]
 def walk(self,tag=None):
  for c in self.children:
   if isinstance(c,Node):
    if tag is None or c.tag==tag:yield c
    yield from c.walk(tag)
 def text(self):return ' '.join(' '.join(c.text() if isinstance(c,Node) else c for c in self.children).split())
class Parser(HTMLParser):
 def __init__(self):super().__init__();self.root=Node();self.stack=[self.root]
 def handle_starttag(self,tag,attrs):
  n=Node(tag,attrs);self.stack[-1].children.append(n)
  if tag not in ['meta','img','br','hr','link','input','source','wbr','area','base','embed','param','col']:self.stack.append(n)
 def handle_startendtag(self,tag,attrs):self.handle_starttag(tag,attrs);self.handle_endtag(tag)
 def handle_endtag(self,tag):
  for i in range(len(self.stack)-1,0,-1):
   if self.stack[i].tag==tag:self.stack=self.stack[:i];break
 def handle_data(self,s):self.stack[-1].children.append(s)
def parse(p, output_dir):
 q=Parser();q.feed(p.read_text());root=q.root
 tables=[]
 for i,t in enumerate(root.walk('table')):
  if 'description-table' in t.attrs.get('class',''):
   rows=[]
   for tr in t.walk('tr'):
    cells=[{'text':x.text(),'images':[im.attrs.get('src') for im in x.walk('img')]} for x in tr.children if isinstance(x,Node) and x.tag in ['td','th']]
    if cells:rows.append(cells)
   tables.append({'table_index':i,'id':t.attrs.get('id'),'text':t.text(),'rows':rows})
 h=p.read_text();family=re.search(r'<h2>ID=(\d+)</h2>',h)
 protocols=[{'id':n.attrs.get('id'),'text':n.text()} for n in root.walk('div') if n.attrs.get('class')=='description-paragraph' and any(v in n.text() for v in ['75 nM','IC50 values','IC 50 values','IC 50 s','IC50s','human PRMT5','Human PRMT5','PRMT5 and MEP50','MEP50 were','PRMT5/MEP50 was'])]
 r={'patent':p.stem,'source':'https://patents.google.com/patent/'+p.stem+'/en','sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'family':family.group(1) if family else None,'meta':{n.attrs.get('name',n.attrs.get('itemprop','')):n.attrs.get('content') for n in root.walk('meta')},'country_representatives':[n.text() for n in root.walk('span') if n.attrs.get('itemprop')=='representativePublication'],'protocols':protocols,'tables':tables}
 (output_dir/(p.stem+'.json')).write_text(json.dumps(r,ensure_ascii=False,indent=2));(output_dir/(p.stem+'.txt')).write_text(root.text());return r
if __name__ == '__main__':
 import argparse
 parser = argparse.ArgumentParser(description="Extract original patent HTML tables with standard-library HTML parsing.")
 parser.add_argument('--raw-dir', type=Path, required=True)
 parser.add_argument('--output-dir', type=Path, required=True)
 args = parser.parse_args()
 args.output_dir.mkdir(parents=True, exist_ok=True)
 for path in sorted(args.raw_dir.glob('*.html')):
  result = parse(path, args.output_dir)
  print(result['patent'], result['family'], len(result['tables']))
