// A completed computation alone is not a passing scientific assessment.
function designRoutePassed(p){
 return !p.running&&!p.is_mock&&p.execution_status==='COMPLETED'&&p.validation==='ADVANCE'&&p.critic==='ADVANCE';
}
function binderLead(p,rows){
 if(!designRoutePassed(p))return null;
 const ranking=(p.tables||[]).find(t=>t.title.startsWith('Additional quality filters'))?.rows||[];
 const rank=new Map(ranking.map(r=>[r.candidate_id,r]));
 const candidates=rows.filter(r=>['PASS','ADVANCE'].includes(r.verdict));
 candidates.sort((a,b)=>{
  const ar=rank.get(a.candidate_id),br=rank.get(b.candidate_id);
  const quality=Number(br?.quality_status==='PASS')-Number(ar?.quality_status==='PASS');
  if(quality)return quality;
  const order=(ar?.rank??Infinity)-(br?.rank??Infinity);
  if(!Number.isNaN(order)&&order!==0)return order;
  return (typeof b.iptm==='number'?b.iptm:-1)-(typeof a.iptm==='number'?a.iptm:-1)||(typeof b.ptm==='number'?b.ptm:-1)-(typeof a.ptm==='number'?a.ptm:-1);
 });
 const row=candidates[0];return row?{row,ranking:rank.get(row.candidate_id)}:null;
}
function renderReportNavigation(reference,current){
 const header=document.querySelector('.topbar');if(!header||document.getElementById('report-navigation'))return;
 const nav=document.createElement('nav');nav.id='report-navigation';nav.className='pr-report-navigation';nav.setAttribute('aria-label','리포트 선택');
 const add=(label,url,active)=>{const a=document.createElement('a');a.textContent=label;a.className='btn ghost';if(url&&/^(?:\/[^/]|\.\.?\/)/.test(url))a.href=url;else a.hidden=true;if(active)a.setAttribute('aria-current','page');nav.append(a);return a;};
 add('표적 발굴 리포트','/report',current==='pipeline');
 const erbb2=add('ERBB2 설계 리포트',reference?.report_url,current==='design');header.prepend(nav);
 if(!window.__COMPETITION_RUN__)fetch('/api/events?cursor=0').then(r=>r.ok?r.json():null).then(d=>{
  if(d?.erbb2_report_run_id){erbb2.href='/report?run_id='+encodeURIComponent(d.erbb2_report_run_id);erbb2.hidden=false;}
 }).catch(()=>{});
}

// Older running servers may not yet include display_notice in figure metadata.
// Read the existing CSV artifact without restarting or rerunning the analysis.
function readReportCSV(text){
 const rows=[];let row=[],cell='',quoted=false;
 for(let i=0;i<text.length;i++){
  const c=text[i];
  if(c==='"'){if(quoted&&text[i+1]==='"'){cell+='"';i++;}else quoted=!quoted;}
  else if(!quoted&&(c===','||c==='\n'||c==='\r')){
   row.push(cell);cell='';if(c!==','){rows.push(row);row=[];if(c==='\r'&&text[i+1]==='\n')i++;}
  }else cell+=c;
 }
 if(cell||row.length){row.push(cell);rows.push(row);}
 const header=rows.shift()||[];return rows.map(values=>Object.fromEntries(header.map((key,i)=>[key,values[i]])));
}
async function explainSequenceFigures(entries,artifacts){
 const sources=new Map();
 for(const [a,figure] of entries){
  const match=String(a.path).match(/^(.*\/)sequence_variation_(\d+)\.png$/i);
  if(!figure||a.display_notice||!match)continue;
  const source=artifacts.find(x=>x.path===match[1]+'sequence_variation.csv'&&x.available&&x.href);
  if(!source)continue;
  if(!sources.has(source.href))sources.set(source.href,[]);
  sources.get(source.href).push([Number(match[2]),figure]);
 }
 await Promise.all([...sources].map(async([href,figures])=>{
  try{
   const response=await fetch(href);if(!response.ok)return;
   const rows=readReportCSV(await response.text());
   for(const [group,figure] of figures){
    const values=rows.filter(r=>Number(r.group)===group);if(!values.length)continue;
    const counts=values.map(r=>Number(r.sequence_count));let notice='';
    if(counts.every(n=>n===1))notice='비교할 서열 부족 · 이 backbone에는 고유 서열이 1개만 있어 서열 간 변화량을 비교할 수 없습니다.';
    else if(counts.every(n=>n>=2)&&values.every(r=>r.mutation_frequency!==''&&Number(r.mutation_frequency)===0))notice='모든 위치의 변화 빈도가 0입니다. 저장된 서열 간 차이가 없습니다.';
    const viewport=figure.querySelector('.report-figure-viewport');
    if(notice&&viewport){const p=document.createElement('p');p.className='ir-figure-notice';p.textContent=notice;viewport.replaceChildren(p);}
   }
  }catch(error){/* Keep the recorded figure when its source table cannot be read. */}
 }));
}

function renderReportFigure(parent,artifact,options={}){
 const el=(host,tag,text,cls)=>{const node=document.createElement(tag);if(text)node.textContent=text;if(cls)node.className=cls;host.append(node);return node;};
 const title=artifact.title||'Analysis Results';
 if(/(^|\/)before_batch_/i.test(artifact.path||''))return null;
 // Pair only matching batch-correction artifacts from the same analysis run.
 const source=String(artifact.path||'');
 const match=source.match(/^(.*\/)?(before|after)_batch_(pca_group|pca_tss|sample_distance)\.[^.]+$/i);
 let host=parent;
 if(match){
  const key=(match[1]||'').replace(/(?:01_qc|02_batch)\/$/i,'')+match[3].toLowerCase();
  host=Array.from(parent.children).find(node=>node.classList.contains('report-comparison')&&node.dataset.comparison===key);
  if(!host){
   host=el(parent,'div',null,'report-comparison');host.dataset.comparison=key;
   host.setAttribute('role','group');host.setAttribute('aria-label','Before and after batch correction');
  }
 }
 const figure=el(host,'figure',null,'report-figure');
 if(match){
  figure.dataset.phase=match[2].toLowerCase();
  if(figure.dataset.phase==='before')host.prepend(figure);
 }

 const header=el(figure,'div',null,'report-figure-header');
 el(header,'h4',title);
 if(artifact.display_notice){el(figure,'p',artifact.display_notice,'ir-figure-notice');return figure;}
 if(artifact.available===false||!artifact.href){el(figure,'p','그림 파일을 찾을 수 없습니다.','ir-figure-notice');return figure;}
 if(options.downloads!==false){const open=el(header,'a','Open full size ↗');
 open.href=artifact.href;open.target='_blank';open.rel='noopener';}
 const viewport=el(figure,'div',null,'report-figure-viewport');
 viewport.tabIndex=0;viewport.setAttribute('role','region');viewport.setAttribute('aria-label',title);
 const link=el(viewport,options.downloads===false?'div':'a',null,'report-figure-link');
 if(options.downloads!==false){link.href=artifact.href;link.target='_blank';link.rel='noopener';link.title='Open full-size image';}
 const img=el(link,'img');img.alt=title;img.loading='lazy';
 img.onerror=()=>{viewport.replaceChildren();el(viewport,'p','그림을 불러오지 못했습니다. 페이지를 새로고침해 다시 시도해 주세요.','ir-figure-notice');};img.src=artifact.href;
 if(artifact.caption)el(figure,'figcaption',artifact.caption);
 return figure;
}

function renderDesignPanel(section,p,options={}){
 const el=(parent,tag,text,cls)=>{const n=document.createElement(tag);if(text!=null)n.textContent=text;if(cls)n.className=cls;parent.append(n);return n;};
 const value=v=>v==null?'—':typeof v==='number'?String(Math.round(v*1000)/1000):typeof v==='object'?JSON.stringify(v):String(v);
 const link=(parent,text,url)=>{const a=el(parent,'a',text);a.href=url;return a;};
 const table=(parent,item)=>{
  el(parent,'h3',item.title);if(!item.rows.length){el(parent,'p','No output recorded for this stage.','ir-muted');return;}
  const wrap=el(parent,'div',null,'ir-scroll'),t=el(wrap,'table'),head=el(t,'tr');item.columns.forEach(key=>el(head,'th',key));
  item.rows.forEach(row=>{const tr=el(t,'tr');if(row===item.lead)tr.className='ir-leading-row';item.columns.forEach(key=>{const cell=el(tr,'td',value(row[key]));if(row===item.lead&&key==='candidate_id')el(cell,'span','대표 후보','ir-result-badge');});});
 };
 const details=(parent,title,data)=>{const d=el(parent,'details');el(d,'summary',title);el(d,'pre',JSON.stringify(data,null,2));return d;};
  el(section,'h2',p.title);if(designRoutePassed(p))el(section,'span','계산 검증 통과','ir-result-badge');el(section,'p',p.scope,'ir-muted');
  el(section,'p','Execution '+p.execution_status+' · Validation '+p.validation+' · Critic '+p.critic);
  el(section,'p','Source '+p.execution_mode+' · Mock '+p.is_mock+' · '+(p.run_id||'Waiting for run'),'ir-muted');
  if(options.downloads!==false&&p.report_url)link(section,'Open modality report ↗',p.report_url);
  if(options.downloads!==false)p.artifacts.filter(a=>a.available&&a.path.endsWith('/binder_results.zip')).forEach(a=>link(section,' · Download analysis ZIP',a.href));
  for(const note of p.analysis_notes||[])el(section,'p',note,'ir-muted');
  const binder=p.modality==='DE_NOVO_BINDER';
  const candidateTables=binder?p.tables.filter(t=>t.title==='Binder candidates'):[];
  const candidateRows=candidateTables.flatMap(t=>t.rows);
  const candidateIds=new Set(candidateRows.map(r=>r.candidate_id));
  // Candidate metrics belong in one comparison table, not repeated cards and bars.
  const summary=p.highlights.filter(m=>!candidateRows.length||(!candidateIds.has(m.label)&&!['iPTM','pTM','Interface contacts'].includes(m.label)));
  if(summary.length){const highlights=el(section,'div',null,'ir-metrics');summary.forEach(m=>{const box=el(highlights,'div');el(box,'small',m.label);el(box,'strong',value(m.value));});}
  const lead=binderLead(p,candidateRows);
  if(lead){
   const callout=el(section,'div',null,'ir-leading-result');el(callout,'strong','대표 Binder 후보 · '+lead.row.candidate_id);
   el(callout,'p','선정 기준: 기본 판정 통과 → 추가 품질 필터 통과 우선 → 저장된 순위 (순위가 없거나 같으면 iPTM·pTM 순).','ir-muted');
   if(lead.ranking)el(callout,'p','추가 품질 필터: '+value(lead.ranking.quality_status)+(lead.ranking.filter_reasons?.length?' · '+lead.ranking.filter_reasons.join(', '):''));
  }
  if(candidateRows.length)table(section,{title:'Binder candidates',rows:lead?[lead.row,...candidateRows.filter(r=>r!==lead.row)]:candidateRows,columns:candidateTables[0].columns,lead:lead?.row});
  table(section,{title:'Step-by-step results',rows:p.stages,columns:['stage','execution_status','message']});
  // Structure visualization is provided by Structure Lab. Reports show recorded metrics and figures.
  const figures=binder?el(section,'div',null,'ir-figures'):section;
  const seenFigures=new Set();
  const uniqueFigures=p.figures.filter(a=>{const key=a.path||a.href;if(!key)return true;if(seenFigures.has(key))return false;seenFigures.add(key);return true;});
  const variations=binder?uniqueFigures.filter(a=>/(^|\/)sequence_variation_\d+\.png$/i.test(a.path||'')):[];
  variations.sort((a,b)=>a.path.localeCompare(b.path,undefined,{numeric:true}));
  const variationSet=new Set(variations);
  const confidenceGroups=new Map();
  uniqueFigures.filter(a=>!variationSet.has(a)).forEach(a=>{
   const match=binder&&String(a.path||'').match(/^(.*\/)binder_(\d+)_(confidence|pae)\.png$/i);
   if(!match){renderReportFigure(figures,a,options);return;}
   const key=match[1]+match[2];
   if(!confidenceGroups.has(key))confidenceGroups.set(key,{label:candidateRows[Number(match[2])]?.candidate_id||'Binder '+match[2],figures:[]});
   confidenceGroups.get(key).figures.push(a);
  });
  if(confidenceGroups.size){
   const detail=el(figures,'details',null,'ir-confidence-figures');el(detail,'summary','AF3 confidence · PAE 그림');
   const label=el(detail,'label','후보 '),select=el(label,'select');select.setAttribute('aria-label','AF3 confidence candidate');
   const groups=[...confidenceGroups.values()];
   groups.forEach((g,i)=>{const option=el(select,'option',g.label);option.value=String(i);});
   const display=el(detail,'div',null,'ir-figures');
   const show=()=>{display.replaceChildren();groups[Number(select.value)].figures.forEach(a=>renderReportFigure(display,a,options));};
   select.onchange=show;show();
  }
  if(variations.length){
   const group=el(figures,'div',null,'ir-sequence-variation');
   el(group,'h3','Binder Sequence Variation');
   const label=el(group,'label','Backbone '),select=el(label,'select');
   select.setAttribute('aria-label','Binder sequence variation backbone');
   const counts=new Map();
   variations.forEach(a=>{const n=a.path.match(/sequence_variation_(\d+)\.png$/i)[1];counts.set(n,(counts.get(n)||0)+1);});
   variations.forEach((a,i)=>{const n=a.path.match(/sequence_variation_(\d+)\.png$/i)[1];const option=el(select,'option','Backbone '+n+(counts.get(n)>1?' · '+a.path.split('/').slice(-3,-1).join('/') : ''));option.value=String(i);});
   const display=el(group,'div',null,'ir-sequence-display');
   const show=()=>{const a=variations[Number(select.value)];display.replaceChildren();const figure=renderReportFigure(display,{...a,title:select.selectedOptions[0].textContent},options);explainSequenceFigures([[a,figure]],p.artifacts);};
   select.onchange=show;show();
  }
  // Keep per-sample provenance accessible without repeating every candidate table.
  const sampleRows=[];let candidateId=null;
  p.tables.forEach(item=>{
   if(binder&&item.title==='Binder candidates'){candidateId=item.rows[0]?.candidate_id;return;}
   if(binder&&item.title==='AF3 samples'){sampleRows.push(...item.rows.map(row=>({...row,candidate_id:row.candidate_id||candidateId})));return;}
   table(section,item);
  });
  if(sampleRows.length){const detail=el(section,'details',null,'ir-af3-samples');el(detail,'summary','AF3 샘플 상세 · '+sampleRows.length+'건');
   const columns=p.tables.find(t=>t.title==='AF3 samples').columns;
   table(detail,{title:'AF3 samples',rows:sampleRows,columns:['candidate_id',...columns.filter(c=>c!=='candidate_id')]});}
  if(p.blockers.length){el(section,'h3','Limitations / reasons');const list=el(section,'ul');p.blockers.forEach(reason=>el(list,'li',reason));}
  if(p.missing_steps.length)el(section,'p','Not evaluated: '+p.missing_steps.join(', '));
  for(const [agent,note] of Object.entries(p.agent_notes||{})){el(section,'h3',agent+' interpretation');el(section,'p',note.text||note.reason||'');}
  details(section,'Execution timeline',p.events);
  if(options.downloads!==false){const files=el(section,'details');el(files,'summary','Files · '+p.artifacts.length);
  p.artifacts.forEach(a=>{const line=el(files,'div',null,'ir-file');if(a.available)link(line,a.label||a.path.split('/').pop(),a.href);else el(line,'span',(a.label||a.path.split('/').pop())+' · unavailable');});}
}

function renderIntegratedRun(data){
 const report=data.integrated_report;
 if(!report)return false;
 const host=document.querySelector('.report-body')||document.querySelector('.content');
 host.replaceChildren();host.classList.add('integrated-report');
 document.title='OmiCraft · ERBB2 Results';
 const title=document.querySelector('.title');if(title)title.textContent='ERBB2 · Therapeutic Design Report';
 const question=document.getElementById('rqText');if(question)question.textContent='Four modalities · Recorded computation and candidate assessment';
 const badge=document.getElementById('verdictBadge');if(badge){badge.textContent=report.workflow_status;badge.className='verdict';}
 const toc=document.querySelector('.toc');if(toc)toc.hidden=true;
 const el=(parent,tag,text,cls)=>{const n=document.createElement(tag);if(text!=null)n.textContent=text;if(cls)n.className=cls;parent.append(n);return n;};
 const link=(parent,text,url)=>{const a=el(parent,'a',text);a.href=url;return a;};
 const top=el(host,'section',null,'ir-hero');
 el(top,'div','ERBB2 / FOUR MODALITIES','ir-eyebrow');el(top,'h2',report.workflow_status==='RUNNING'?'Live results':'Design results');
 el(top,'p',`${report.completed_routes} / 4 routes ended · Workflow ${report.workflow_status} · Report ${report.report_status}`);
 const progress=el(top,'progress');progress.max=4;progress.value=report.completed_routes;
 el(top,'p',data.run_id+' · '+(data.execution_mode||'Recorded results'),'ir-muted');
 el(top,'p','계산의 종료 상태와 후보의 과학적 판정은 별도로 표시합니다. 실패하거나 미평가인 단계도 아래에 남습니다.','ir-muted');
 const actions=el(top,'div',null,'ir-actions');link(actions,'Earlier ERBB2 Report ↗',report.previous_report_url);
 const refresh=el(actions,'button','Refresh');refresh.onclick=()=>location.reload();
 if(report.passing_candidates.length){
  const callout=el(host,'section',null,'ir-success');el(callout,'h3','Structural screening passed');
  report.passing_candidates.forEach(p=>{const a=link(callout,p.modality.replaceAll('_',' ')+' · ADVANCE ↗',p.report_url);a.style.display='block';});
  el(callout,'p','현재 계산 기준을 통과한 후보입니다. 실험적 결합·약효를 의미하지 않습니다.');
 }
 const cards=el(host,'div',null,'ir-cards');const navigation=el(host,'nav',null,'ir-tabs');
 const sections=[];
 function select(name){sections.forEach(([key,node])=>node.hidden=name!=='all'&&key!==name);navigation.querySelectorAll('button').forEach(b=>b.classList.toggle('selected',b.dataset.tab===name));window.omicraftSelectedTab=name;}
 const all=el(navigation,'button','All results');all.dataset.tab='all';all.onclick=()=>select('all');
 for(const p of report.panels){
  const card=el(cards,'button',null,'ir-card');card.onclick=()=>select(p.modality);el(card,'strong',p.title);el(card,'span',p.execution_status,'ir-status '+(p.validation==='ADVANCE'?'pass':''));el(card,'small','Assessment · '+p.validation);
  const tab=el(navigation,'button',p.title);tab.dataset.tab=p.modality;tab.onclick=()=>select(p.modality);
  const section=el(host,'section',null,'ir-section');section.id='result-'+p.modality;sections.push([p.modality,section]);
  renderDesignPanel(section,p);
 }
 select(window.omicraftSelectedTab||'all');
 if(report.workflow_status==='RUNNING'&&!window.__COMPETITION_RUN__){
  clearTimeout(window.omicraftReportTimer);window.omicraftReportTimer=setTimeout(async()=>{
   try{const response=await fetch('/api/report'+location.search,{cache:'no-store'});if(response.ok)window.renderTherapeuticRun(await response.json());}catch(error){console.warn('Report refresh:',error.message);}
  },10000);
 }
 return true;
}
