// Optional Node.js DOM-harness smoke: saved choices, traces and nearest-row hover.
// This tests control logic, not browser rendering or accessibility.
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import assert from 'node:assert/strict';
import {fileURLToPath} from 'node:url';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const html=fs.readFileSync(path.join(root,'docs/model_viewer.html'),'utf8');
const raw=html.match(/<script id="saved-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const code=html.match(/<script>\s*([\s\S]*?)<\/script>/)[1];
const saved=JSON.parse(raw);
class Element {
  constructor(){this.children=[];this.events={};this.attributes={};this.style={};this.value='';this.checked=false;this.textContent='';}
  add(option){this.children.push(option);if(this.value==='')this.value=String(option.value);}
  append(...children){this.children.push(...children);}
  appendChild(child){this.children.push(child);return child;}
  replaceChildren(...children){this.children=[...children];this.value='';}
  setAttribute(name,value){this.attributes[name]=value;}
  addEventListener(name,callback){this.events[name]=callback;}
  getBoundingClientRect(){return {left:0,width:1000};}
}
const ids=Object.fromEntries(['saved-data','model','policy','compare','chart','equation','context','legend','source','sample'].map(id=>[id,new Element()]));
ids['saved-data'].textContent=raw;
const document={getElementById:id=>ids[id],createElement:()=>new Element(),createElementNS:()=>new Element(),createTextNode:text=>({textContent:text})};
function Option(text,value){this.text=text;this.value=String(value);}
vm.runInNewContext(code,{document,Option,console},{timeout:5000});
let choices=0,comparisons=0,hovers=0;
for(let i=0;i<saved.length;i++){
  ids.model.value=String(i);ids.model.events.change();
  assert.equal(ids.policy.children.length,Object.keys(saved[i].series).length);
  assert.ok(ids.equation.textContent.includes(saved[i].model));
  assert.ok(ids.source.textContent.includes(saved[i].sha256));
  for(const key of Object.keys(saved[i].series)){
    ids.policy.value=key;ids.policy.events.change();
    const traces=ids.chart.children.filter(n=>n.attributes.points);
    assert.equal(traces.length,['M3','M4'].includes(saved[i].model)?2:1);
    assert.ok(traces.every(n=>!n.attributes.points.includes('NaN')&&!n.attributes.points.includes('null')));
    ids.chart.events.pointermove({clientX:520});
    assert.ok(ids.sample.textContent.startsWith(key+' | t='));
    assert.ok(ids.sample.textContent.includes(' | x='));
    assert.equal(ids.sample.textContent.includes(' | s='),saved[i].model==='M3');
    assert.equal(ids.sample.textContent.includes(' | d='),saved[i].model==='M4');
    choices++;hovers++;
  }
  ids.compare.checked=true;ids.compare.events.change();
  assert.equal(ids.chart.children.filter(n=>n.attributes.points).length,Object.keys(saved[i].series).length);
  comparisons++;
}
assert.equal(choices,37);
console.log(JSON.stringify({datasets:saved.length,policy_choices:choices,comparisons,nearest_row_hover_checks:hovers,status:'passed',scope:'DOM logic only; browser rendering not tested'}));
