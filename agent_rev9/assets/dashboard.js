const $ = id => document.getElementById(id);
const icon = name => `<svg class="icon" aria-hidden="true"><use href="#i-${name}"/></svg>`;
const stages = [
 ['planner','Planner','Research question & contrasts'],
 ['discovery','Discovery','Omics evidence & interpretation'],
 ['qualification','Qualification','Tier assessment & safety gates'],
 ['design','Design','Modality-specific workflows'],
 ['critic','Critic','Evidence review'],
 ['dossier','Dossier','Results & next steps']
];
const groups = {resolve_contrasts_node:'planner',planner_review:'planner',planner:'planner',r_analysis_node:'discovery',apear_node:'discovery',cell_context_node:'discovery',depmap_node:'discovery',discovery_review:'discovery',discovery:'discovery',r_analysis:'discovery',apear:'discovery',cell_context:'discovery',depmap:'discovery',qualification_node:'qualification',qualification_review:'qualification',assess_modalities:'qualification',user_selection_gate_node:'qualification',user_gate:'qualification',automatic_design_node:'design',design:'design',therapeutic_design:'design',critic_node:'critic',critic_review:'critic',critic:'critic',output_tables_node:'qualification',output_tables:'qualification',dossier_review:'dossier'};
const terminalNodes = {planner_review:'planner',planner:'planner',discovery_review:'discovery',discovery:'discovery',qualification_review:'qualification',qualification:'qualification',automatic_design_node:'design',design:'design',critic_review:'critic',critic:'critic',dossier_review:'dossier'};
const nodeAlias = {planner:'resolve_contrasts_node',r_analysis:'r_analysis_node',apear:'apear_node',cell_context:'cell_context_node',depmap:'depmap_node',qualification:'qualification_node',output_tables:'output_tables_node',critic:'critic_review'};
const state = {runId:0,cursor:0,running:false,connected:false,requesting:false,polling:false,report:'',started:0,events:[],logs:[],filter:'all',targets:[],choices:[],catalog:[],catalogReady:false,allTargets:false,steps:{},questionEdited:false,qualification:false};
const emptySteps = () => Object.fromEntries(stages.map(([name])=>[name,{status:'waiting',elapsed:null}]));
state.steps=emptySteps();
function element(tag,text,cls){const node=document.createElement(tag);if(text!==undefined)node.textContent=text;if(cls)node.className=cls;return node;}
function clean(value){return String(value??'').replace(/\/(?:NHNHOME|home|tmp)\/[^\s"',;]+/g,path=>path.split('/').pop());}
function groupFor(node){return node.startsWith('R:')?'discovery':/^(small_molecule|de_novo_binder|adc|degrader)\//.test(node)?'design':groups[node]||node;}
function controls(){const disabled=!state.connected||state.running||state.requesting;$('btnRun').disabled=disabled;$('btnERBB2').disabled=disabled;$('btnStop').hidden=!state.running;$('btnStop').disabled=!state.connected||state.stopping||state.cancelling;$('btnStop').textContent=state.cancelling?'중단 중…':'Stop Pipeline';}
function renderWorkflow(){
 $('workflow').replaceChildren();let completed=0;
 for(const [name,label,description] of stages){
  const step=state.steps[name];if(step.status==='done')completed++;
  const item=element('li',undefined,step.status==='done'?'done':step.status==='running'?'active-step':step.status==='failed'?'step-failed':'');
  const circle=element('span',undefined,'step-circle');
  if(step.status==='done')circle.innerHTML=icon('check');else circle.textContent=String(stages.findIndex(s=>s[0]===name)+1).padStart(2,'0');
  const copy=element('span',undefined,'step-copy');copy.append(element('span',label,'step-name'),element('p',description,'step-desc'));
  item.append(circle,copy,element('span',step.elapsed==null?'':step.elapsed.toFixed(1)+'s','step-time mono'));
  $('workflow').append(item);
 }
 $('stepsLabel').textContent=completed+' / 6 stages';
 renderDetailedProgress();
}
const detailStages = [
 ['resolve_contrasts_node','Contrasts','Contrast specification'],
 ['r_analysis_node','R analysis','DESeq2 · GSEA'],
 ['apear_node','GSEA','Pathway enrichment'],
 ['cell_context_node','Cell context','Cell-of-origin evidence'],
 ['depmap_node','DepMap','Dependency evidence'],
 ['qualification_node','Qualification','Target qualification'],
 ['output_tables_node','Tables','Evidence tables'],
 ['assess_modalities','Modalities','Modality assessment'],
 ['user_selection_gate_node','Gate','Design eligibility'],
 ['automatic_design_node','Design','Therapeutic design'],
 ['critic_review','Critic','Evidence review'],
 ['dossier_review','Dossier','Results & next steps']
];
let detailState, detailCurrent, detailTerminal, detailRoutes;
function resetDetailedProgress(){
 detailState=Object.fromEntries(detailStages.map(([name])=>[name,{status:'waiting',elapsed:null,duration:null,start:null,message:''}]));
 detailCurrent=null;detailTerminal=null;detailRoutes=new Map();
}
function visibleDetailStages(){
 if(state.report.startsWith('lab_'))return detailStages.filter(([name])=>name==='automatic_design_node');
 if(state.report.startsWith('design_'))return detailStages.filter(([name])=>['automatic_design_node','critic_review'].includes(name));
 return detailStages;
}
function detailStatus(event){
 const status=String(event.execution_status||'');
 if(['HOLD','REJECT'].includes(status))return 'blocked';
 if(['PENDING','UNKNOWN'].includes(status))return 'pending';
 if(status==='CANCELLED')return 'cancelled';
 if(status.startsWith('BLOCKED'))return 'blocked';
 if(status==='FAILED'||(event.errors||[]).length)return 'failed';
 if(['NOT_RUN','NOT_REQUESTED','SKIPPED'].includes(status))return 'skipped';
 if(['PARTIAL','COMPLETED_WITH_ERRORS'].includes(status))return 'partial';
 return 'done';
}
function detailedEvent(event){
 const node=String(event.node||'');
 // LangGraph combines modality assessment and eligibility in one user_gate node.
 if(node==='user_gate'&&(event.type==='node'||event.type==='progress')){
  detailedEvent({...event,node:'assess_modalities'});
  detailedEvent({...event,node:'user_selection_gate_node',execution_status:event.gate_status||'UNKNOWN'});
  return;
 }
 const isR=node.startsWith('R:');
 const isRoute=/^(small_molecule|de_novo_binder|adc|degrader)\//.test(node);
 const isDesign=node==='therapeutic_design'||node==='design'||isRoute;
 const name=isR?'r_analysis_node':isDesign?'automatic_design_node':node==='critic_node'&&state.report.startsWith('design_')?'critic_review':(nodeAlias[node]||node);
 const item=detailState[name];
 const messages=[...(event.messages||[]),...(event.errors||[])].filter(Boolean);
 const message=clean(messages.at(-1)||'').replace(/^\[[^\]]+\]\s*/,'');
 if(event.type==='progress'||event.type==='node'){
  if(item){
   if(event.type==='progress'){
    if(item.status!=='running'){item.start=isDesign?null:(event.elapsed??null);item.duration=null;}
    item.status='running';
   }else if(!isR&&!isRoute&&node!=='therapeutic_design'){
    item.status=detailStatus(event);item.elapsed=event.elapsed??null;
    if(item.start!=null&&event.elapsed>=item.start)item.duration=event.elapsed-item.start;
   }
   if(isR)item.message=node.slice(2).trim();else if(message)item.message=message;
   detailCurrent={name,label:detailStages.find(([key])=>key===name)?.[2]||node,message:item.message,status:event.execution_status||item.status};
   if(isDesign&&message){
    const routeName=message.match(/^([A-Za-z0-9_.-]+:(?:SMALL_MOLECULE|DE_NOVO_BINDER|ADC|DEGRADER))\b/)?.[1];
    const key=routeName||isRoute&&node;
    if(key)detailRoutes.set(key,{label:routeName||node.replaceAll('_',' ').replace('/',' · '),message,status:event.type==='progress'?'running':detailStatus(event)});
   }
  }else if(['planner_review','discovery_review','qualification_review'].includes(node)){
   detailCurrent={name:null,label:node.replace('_review','')+' · interpretation',message,status:event.type==='progress'?'RUNNING':'COMPLETED'};
  }
 }else if(event.type==='done'||event.type==='error'){
  detailTerminal=event;
  const terminalStatus=event.type==='error'?'failed':detailStatus(event);
  for(const step of Object.values(detailState)){
   if(step.status==='running')step.status=terminalStatus;
   else if(step.status==='waiting')step.status='skipped';
  }
  for(const route of detailRoutes.values())if(route.status==='running')route.status=terminalStatus==='done'?'unknown':terminalStatus;
 }
}
function renderDetailedProgress(){
 const rows=visibleDetailStages();
 const finished=rows.filter(([name])=>['done','partial'].includes(detailState[name].status)).length;
 const skipped=rows.filter(([name])=>detailState[name].status==='skipped').length;
 const percent=Math.round(finished/rows.length*100);
 const running=state.running&&!detailTerminal;
 $('detailCount').textContent=finished+' / '+rows.length+' steps'+(skipped?' · '+skipped+' not run':'');
 $('detailPercent').textContent=percent+'%';
 $('progressFill').style.width=percent+'%';
 $('progressFill').className='progress-fill'+(running?' is-running':'')+(detailTerminal?.type==='error'?' is-failed':'');
 $('stageProgress').setAttribute('aria-valuenow',percent);
 $('stageProgress').setAttribute('aria-valuetext',finished+' of '+rows.length+' steps completed');
 $('detailSegments').replaceChildren();
 for(const [name,label,description] of rows){
  const item=detailState[name];const index=detailStages.findIndex(([key])=>key===name)+1;
  const segment=element('div',undefined,'detail-segment '+item.status);
  segment.setAttribute('title','WF'+String(index).padStart(2,'0')+' · '+description+' · '+item.status+(item.message?' — '+item.message:''));
  if(item.status==='running')segment.setAttribute('aria-current','step');
  segment.append(element('span',item.status==='done'?'✓':String(index).padStart(2,'0'),'detail-circle'),element('span',label,'detail-name'));
  const elapsed=item.duration!=null?item.duration.toFixed(1)+'s':item.elapsed!=null?'@ '+Number(item.elapsed).toFixed(1)+'s':item.status==='skipped'?'Not run':item.status==='running'?'Running':item.status==='waiting'?'':item.status;
  segment.append(element('span',elapsed,'detail-time mono'));$('detailSegments').append(segment);
 }
 const current=detailCurrent;
 $('detailActivity').className='detail-activity'+(running?' is-running':'');
 $('detailCurrent').textContent=detailTerminal?String(detailTerminal.execution_status||'FAILED'):current?current.label:'Ready';
 $('detailMessage').textContent=detailTerminal?'':current?.message||'';
 $('detailRoutes').replaceChildren();
 for(const route of detailRoutes.values()){
  const row=element('div',undefined,'detail-route '+route.status);row.setAttribute('title',route.message);
  row.append(element('span',route.label),element('span',route.status==='done'?'Completed':route.status==='unknown'?'Status not reported':route.status.replace(/^./,c=>c.toUpperCase()),'mono'));
  $('detailRoutes').append(row);
 }
 $('detailRoutes').hidden=detailRoutes.size===0;
}
resetDetailedProgress();

function appendLog(event,message,group,error=false){
 if(!String(message??'').trim())return;
 const text=clean(message);const match=text.match(/^\[([^\]]+)\]\s*/);
 state.logs.push({time:event.elapsed??event.total_time??null,group,label:match?match[1]:(stages.find(s=>s[0]===group)?.[1]||group||'Pipeline'),text:match?text.slice(match[0].length):text,error});
}
function renderLog(){
 const visible=state.logs.filter(line=>state.filter==='all'||line.group===state.filter);
 $('eventCount').textContent=visible.length+' events';$('logWindow').replaceChildren();
 for(const entry of visible.slice(-2000)){
  const row=element('div',undefined,'log-line '+entry.group+(entry.error?' error':''));
  row.append(element('span',entry.time==null?'':Number(entry.time).toFixed(1)+'s','ts'),element('b','['+entry.label+']'),element('span',entry.text,'log-text'));
  $('logWindow').append(row);
 }
 if(!visible.length)$('logWindow').append(element('p','No events for this stage.','small-note'));
 if($('autoScroll').checked)$('logWindow').scrollTop=$('logWindow').scrollHeight;
}
function updateReportEntry(){
 const kind=$('reportKind');if(!kind)return;
 const erbb2=kind.querySelector('[value="erbb2"]');
 erbb2.disabled=!state.erbb2Report;erbb2.textContent=state.erbb2Report?'ERBB2 · 모달리티 설계':'ERBB2 · 설계 결과 불러오는 중';
 $('pipelineReport').href=kind.value==='erbb2'&&state.erbb2Report?'/report?run_id='+encodeURIComponent(state.erbb2Report):'/report';
}
if($('reportKind'))$('reportKind').onchange=updateReportEntry;
function reportLink(run){
 if(!run)return;
 const current=$('currentReport');current.href='/report?run_id='+encodeURIComponent(run);current.hidden=false;
}
async function loadERBB2(run){
 try{
  const url=run?'/report?run_id='+encodeURIComponent(run):$('savedERBB2').getAttribute('href');
  const report=await api(url.replace('/report?','/api/report?'));
  if(report.gene!=='ERBB2'||report.modality!=='ALL'||(run&&state.erbb2Report&&state.erbb2Report!==run))return;
  $('savedERBB2').href=url;state.erbb2Report=report.run_id;updateReportEntry();
  const integrated=report.integrated_report;
  if(integrated)$('erbb2Summary').textContent=integrated.completed_routes+' / 4 routes ended · '+integrated.passing_candidates.length+' modality passed · '+integrated.workflow_status;
  const cards=$('erbb2Results');cards.replaceChildren();
  for(const row of report.integrated_report?.panels||report.summary?.modalities||[]){
   const child=/^design_\d{8}T\d{6}Z_[0-9a-f]{8}$/.test(row.run_id||'')?'/report?run_id='+row.run_id:url;
   const a=link(row.modality.replaceAll('_',' '),child);
   if((row.validation||row.validation_status)==='ADVANCE')a.classList.add('passed-result');
   a.append(element('small',(row.running?'RUNNING':row.execution_status||'NOT_RUN')+' · '+(row.validation||row.validation_status||row.critic_decision||'NOT_EVALUATED')));cards.append(a);
  }
 }catch(error){console.warn('ERBB2 results:',error.message);}
}
async function loadHistory(){
 try{
  const response=await fetch('/api/reports',{cache:'no-store'});
  if(!response.ok)return;
  const data=await response.json(),select=$('historySelect'),selected=select.value;
  select.replaceChildren(element('option','Select a run'));select.firstChild.value='';
  for(const row of data.runs){
   const option=element('option',[row.run_id.replace(/_[0-9a-f]{8}$/,''),row.gene||'',row.modality||'',row.execution_status].filter(Boolean).join(' · '));
   option.value=row.report_url;select.append(option);
  }
  select.value=selected;select.onchange();$('runHistory').hidden=false;
 }catch(error){console.warn('Local report history:',error.message);}
}
function link(label,href,cls){const node=element('a',label,cls);node.href=href;node.target='_blank';node.rel='noopener';return node;}
function passedDesign(row){return !row.running&&row.execution_status==='COMPLETED'&&row.validation==='ADVANCE'&&row.critic==='ADVANCE'&&!row.is_mock;}
function targetsFromReport(data){
 const hidden=new Set(data.candidate_display?.hidden_genes||[]);
 const panels=data.design_report?.panels||[];
 return (data.summary?.candidates||[]).filter(t=>!hidden.has(t.gene)).map(t=>{
  const results=panels.filter(p=>p.gene===t.gene);
  return {...t,parent:data.run_id,results,designPassed:results.some(passedDesign)};
 });
}
function renderTargets(){
 const final=!!state.finalTargets;
 const targets=state.targets.filter(row=>state.allTargets||(final?row.designPassed:row.decision==='ADVANCE'));
 $('targetList').replaceChildren();$('targetCount').textContent=targets.length+' targets';
 $('targetsTitle').textContent=state.allTargets?'Evaluated Targets':final?'ADVANCE Targets · 설계 검증 통과':'ADVANCE Targets · 설계 전 적격성';
 if(!targets.length)$('targetList').append(element('p',final?'설계 검증을 통과한 표적이 없습니다.':'No targets to display.','small-note'));
 if(state.targetsReport)$('targetList').append(element('p','기준 실행: '+state.targetsReport,'small-note'));
 for(const target of targets){
  const gene=String(target.gene||'').toUpperCase();const card=element('article',undefined,'target-item');
  const header=element('div',undefined,'flex spread');header.append(element('h3',gene,'gene'),element('span',target.tier?'Tier '+target.tier:target.decision,'pill navy'));card.append(header);
  card.append(element('p','설계 전 적격성: '+target.decision,'target-modality'));
  for(const row of target.results||[]){
   if(!state.allTargets&&!passedDesign(row))continue;
   card.append(element('p',row.title+' · '+row.execution_status+' · Validation '+row.validation+' / Critic '+row.critic,'target-modality'));
   if(row.report_url)card.append(link('결과 보고서',row.report_url,'design-link'));
   if((row.structures||[]).some(a=>a.available))card.append(link('Structure Lab','/lab?'+new URLSearchParams({run_id:target.parent,gene,modality:row.modality}),'design-link'));
  }
  $('targetList').append(card);
 }
 $('showAdvance').classList.toggle('active',!state.allTargets);$('showAllTargets').classList.toggle('active',state.allTargets);
}
let targetRequest=0;
async function loadLatestTargets(id){
 if(!id)return;
 const ticket=++targetRequest,runId=state.runId;
 try{
  const data=await api('/api/report?run_id='+encodeURIComponent(id));
  if(ticket!==targetRequest||runId!==state.runId||state.running||data.running)return;
  state.targets=targetsFromReport(data);state.targetsReport=id;state.finalTargets=(data.design_report?.panels||[]).length>0;
  $('metricEvaluated').textContent=state.targets.length;
  $('metricAdvance').textContent=state.targets.filter(t=>state.finalTargets?t.designPassed:t.decision==='ADVANCE').length;
  $('advanceFoot').textContent=state.finalTargets?'validated design targets':'qualified for design';
  $('metricHold').textContent=state.targets.filter(t=>t.decision==='HOLD').length;
  if(data.summary?.dossier?.status==='COMPLETED' && state.report===id){
   state.steps.dossier={status:'done',elapsed:null};detailState.dossier_review={...detailState.dossier_review,status:'done',message:'저장된 Dossier 생성 완료'};renderWorkflow();
  }
  renderTargets();
 }catch(error){console.warn('Latest pipeline results:',error.message);}
}

function resetRun(){
 resetDetailedProgress();
 targetRequest++;state.finalTargets=false;state.targetsReport='';state.targetsRefreshed=0;
 $('advanceFoot').textContent='qualified for design';
 state.logs=[];state.events=[];state.targets=[];state.steps=emptySteps();state.qualification=false;
 for(const id of ['metricEvaluated','metricAdvance','metricHold','scientificVerdict'])$(id).textContent='—';
 $('criticIterations').textContent='0';$('overviewSubtitle').textContent='';$('progressCaption').textContent='';
 renderWorkflow();renderLog();renderTargets();
}
function consume(event){
 detailedEvent(event);
 state.events.push(event);const node=String(event.node||'');const group=groupFor(node);
 const step=state.steps[group];
 if(event.type==='progress'){
  if(step&&node!=='critic_node'){step.status='running';}
  for(const message of event.messages||[])if(node.startsWith('R:')||node==='therapeutic_design'||/^(small_molecule|de_novo_binder|adc|degrader)\//.test(node))appendLog(event,message,group);
  if(step&&node!=='critic_node'){$('overviewTitle').textContent=stages.find(s=>s[0]===group)[1];$('progressCaption').textContent=node.replaceAll('_',' ');}
 }else if(event.type==='node'){
  for(const message of event.messages||[])appendLog(event,message,group);
  for(const issue of event.review_issues||[])appendLog(event,'검토 의견: '+issue,'critic');
  for(const error of event.errors||[])appendLog(event,error,group,true);
  if(terminalNodes[node]){const item=state.steps[terminalNodes[node]];item.status=detailStatus(event);item.elapsed=Number(event.elapsed)||0;}
  if(event.has_qualification_counts){
   state.qualification=true;
   $('metricAdvance').textContent=event.n_advance??0;$('metricHold').textContent=event.n_hold??0;
   $('metricEvaluated').textContent=(event.n_advance||0)+(event.n_hold||0)+(event.n_reject||0);
   state.targets=(event.evaluated_targets||[]).map(t=>({...t,parent:state.report}));
   $('overviewSubtitle').textContent=(event.advance_targets||[]).map(t=>t.gene).join(' · ');
   renderTargets();loadChoices();
  }
  if(event.verdict){$('scientificVerdict').textContent=event.verdict;}
 }else if(event.type==='done'||event.type==='error'){
  state.running=false;
  const cancelled=event.execution_status==='CANCELLED';
  const failed=event.type==='error'||event.execution_status==='FAILED';
  $('runStatus').textContent=cancelled?'중단됨':event.execution_status||(failed?'FAILED':'COMPLETED');
  $('runStatus').className='pill '+(failed||event.execution_status!=='COMPLETED'?'amber':'teal');
  $('scientificVerdict').textContent=event.final_verdict||event.scientific_validation_status||'NOT_EVALUATED';
  $('overviewTitle').textContent=cancelled?'실험 중단됨':'Pipeline results';$('progressCaption').textContent='';
  $('criticIterations').textContent=event.critic_count??0;
  if(event.design_outcome&&Object.keys(event.design_outcome).length)$('overviewSubtitle').textContent='계산 완료 '+event.design_outcome.completed_routes+'개 · 검증 통과 '+event.design_outcome.validated_successes+'개';
  if(event.total_time!=null)$('elapsed').textContent=formatTime(event.total_time);
  for(const error of [...new Set([...(event.errors||[]),...(event.blockers||[])])])appendLog(event,error,'critic',true);
  for(const step of Object.values(state.steps))if(step.status==='running')step.status=detailStatus(event);
  reportLink(event.report_run_id);loadHistory();if(event.erbb2_batch)loadERBB2(event.report_run_id);loadChoices();loadStructures();
 }
}
function formatTime(seconds){const value=Math.max(0,Math.floor(seconds));return String(Math.floor(value/60)).padStart(2,'0')+':'+String(value%60).padStart(2,'0');}
async function api(path,body){
 const response=await fetch(path,{method:body===undefined?'GET':'POST',cache:'no-store',...(body===undefined?{}:{headers:{'Content-Type':'application/json'},body:JSON.stringify(body)})});
 let data;try{data=await response.json();}catch(error){throw new Error('HTTP '+response.status);}
 if(!response.ok)throw new Error(typeof data.detail==='string'?data.detail:'HTTP '+response.status);
 return data;
}
async function loadChoices(){
 try{
  state.choices=(await api('/api/design/candidates')).filter(row=>row.origin==='upstream_analysis');

  renderTargets();
 }catch(error){console.warn('Design candidates:',error.message);}
}
async function loadStructures(){try{state.catalog=await api('/api/structures');state.catalogReady=true;renderTargets();}catch(error){console.warn('Structures:',error.message);}}
async function resources(){
 try{const data=await api('/api/resources');$('gpuRow').replaceChildren();for(const gpu of data.devices){const card=element('div',undefined,'gpu '+(gpu.available===false?'busy':''));card.append(element('strong','GPU '+gpu.index,'mono'),element('span',gpu.available==null?'—':gpu.available?'Available':'In use'));$('gpuRow').append(card);}}
 catch(error){console.warn('Resources:',error.message);}finally{setTimeout(resources,15000);}
}
let pollTimer;
async function poll(){
 if(state.polling)return;clearTimeout(pollTimer);state.polling=true;
 try{
  const data=await api('/api/events?run_id='+state.runId+'&cursor='+state.cursor);
  if(data.run_id!==state.runId||data.cursor<state.cursor){state.runId=data.run_id;state.cursor=0;resetRun();}
  state.connected=true;state.cancelling=!!data.cancelling;state.running=data.running;state.started=data.started_at||0;state.report=data.report_run_id||state.report;
  reportLink(state.report);
  // Pipeline Report uses /report to resolve the latest usable analysis at click time.
  if(data.erbb2_report_run_id&&(data.erbb2_report_run_id!==state.erbb2Report||Date.now()-(state.erbb2Refreshed||0)>5000)){
   state.erbb2Refreshed=Date.now();
   state.erbb2Report=data.erbb2_report_run_id;loadERBB2(state.erbb2Report);
  }
  if(!state.questionEdited&&typeof data.research_question==='string'){
   const option=[...$('researchQuestion').options].find(item=>item.dataset.question===data.research_question);
   if(option){$('researchQuestion').value=option.value;updateQuestion();}
  }
  for(const event of data.events||[])consume(event);
  state.cursor=data.cursor;
  const latest=data.pipeline_report_run_id;
  if(!data.running&&latest&&(latest!==state.targetsReport||Date.now()-(state.targetsRefreshed||0)>10000)){
   state.targetsRefreshed=Date.now();loadLatestTargets(latest);
  }
  if(data.running){
   $('runStatus').textContent=data.cancelling?'중단 중…':'RUNNING';$('runStatus').className='pill teal';
   if(['lab','design'].includes(data.execution_kind)&&data.task_label)$('overviewTitle').textContent=data.task_label;
  }
  else if(!data.run_id){$('runStatus').textContent='Ready';$('overviewTitle').textContent='Research workspace';}
  renderWorkflow();if((data.events||[]).length)renderLog();
 }catch(error){state.connected=false;console.warn('Status connection:',error.message);}
 finally{state.polling=false;controls();pollTimer=setTimeout(poll,1500);}
}
async function start(path,body){
 state.requesting=true;controls();
 try{const result=await api(path,body);if(!['started','already_running'].includes(result.status))throw new Error('Execution was not started');state.running=true;await poll();}
 catch(error){appendLog({},error.message,'planner',true);renderLog();await poll();}
 finally{state.requesting=false;controls();}
}
$('historySelect').onchange=()=>{const url=$('historySelect').value;const button=$('historyReport');button.hidden=!url;if(url)button.href=url;};
loadHistory();
$('btnStop').onclick=async()=>{state.stopping=true;controls();try{await api('/api/stop',{run_id:state.runId});await poll();}catch(error){appendLog({},error.message,'planner',true);renderLog();}finally{state.stopping=false;controls();}};
function showQuestionError(msg){const el=$('questionError');if(!el)return;if(msg){el.textContent=msg;el.style.display='block';}else{el.textContent='';el.style.display='none';}}
$('btnRun').onclick=async()=>{
  const scenario=$('researchQuestion').value;
  if(!scenario){$('researchQuestion').focus();return;}
  showQuestionError('');
  const structural_backend=($('structuralBackend')||{}).value||'af3';
  state.requesting=true;controls();
  try{
    const result=await api('/api/run',{scenario,design_mode:'all_eligible',structural_backend});
    if(!['started','already_running'].includes(result.status))throw new Error('Execution was not started');
    state.running=true;await poll();
  }catch(error){
    showQuestionError(error.message);
    await poll();
  }finally{state.requesting=false;controls();}
};
$('btnERBB2').onclick=()=>start('/api/design',{parent_run_id:'reference_erbb2',gene:'ERBB2',modality:'ALL',input:{}});
function updateQuestion(){
 const select=$('researchQuestion'),option=select.selectedOptions[0];
 $('questionDescription').textContent=option.dataset.question;
 const erbb2=select.value==='erbb2_design';
 $('questionCohort').textContent=erbb2?'ERBB2 · 표적 구조 및 모달리티별 입력':'TCGA-BRCA · 1,018 samples';
 $('questionFlow').textContent=erbb2?'저분자 · 단백질 바인더 · ADC · 단백질 분해제':'Discovery → Qualification → Design';
 $('structuralBackend').disabled=erbb2;
 if(erbb2)$('structuralBackend').value='af3';
}
$('researchQuestion').addEventListener('change',()=>{state.questionEdited=true;showQuestionError('');updateQuestion();});
updateQuestion();
$('showAdvance').onclick=()=>{state.allTargets=false;renderTargets();};$('showAllTargets').onclick=()=>{state.allTargets=true;renderTargets();};
document.querySelectorAll('[data-filter]').forEach(button=>button.onclick=()=>{state.filter=button.dataset.filter;document.querySelectorAll('[data-filter]').forEach(tab=>{const active=tab===button;tab.classList.toggle('selected',active);tab.setAttribute('aria-pressed',String(active));});renderLog();});
$('copyLog').onclick=()=>{const text=state.logs.map(e=>`${e.time??''}s [${e.label}] ${e.text}`).join('\n');const url=URL.createObjectURL(new Blob([text],{type:'text/plain;charset=utf-8'}));const anchor=document.createElement('a');anchor.href=url;anchor.download='omicraft-execution-log.txt';anchor.click();setTimeout(()=>URL.revokeObjectURL(url),1000);};
setInterval(()=>{if(state.running&&state.started)$('elapsed').textContent=formatTime(Date.now()/1000-state.started);},1000);
renderWorkflow();renderLog();renderTargets();controls();poll();loadChoices();loadStructures();resources();
