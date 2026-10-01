(()=>{
'use strict';
const $=id=>document.getElementById(id),ns='http://www.w3.org/2000/svg';
let topology=null,data=null,selected=null,transport=null,connectionEpoch=0,lastKey='',lastRun='',requesting=false;
const text=(tag,value,cls)=>{const n=document.createElement(tag);n.textContent=value??'';if(cls)n.className=cls;return n;};
const svg=(tag,attrs,parent)=>{const n=document.createElementNS(ns,tag);for(const[k,v]of Object.entries(attrs))n.setAttribute(k,v);parent.append(n);return n;};
const duration=ms=>ms==null?'—':`${(ms/1000).toFixed(1)}s`;
const time=value=>value?new Date(value).toLocaleTimeString():'—';
async function api(path,body){const r=await fetch(path,{cache:'no-store',...(body?{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify(body)}:{})});const d=await r.json();if(!r.ok)throw Error(d.detail||`HTTP ${r.status}`);return d;}
function connection(message,error=false){$('connection').textContent=message;$('connection').className='connection'+(error?' error':'');}
function draw(){
 const graph=$('graph');graph.replaceChildren();if(!topology?.nodes?.length){svg('text',{x:30,y:50},graph).textContent='Graph definition unavailable';return;}
 const defs=svg('defs',{},graph),marker=svg('marker',{id:'arrow',viewBox:'0 0 10 10',refX:9,refY:5,markerWidth:6,markerHeight:6,orient:'auto-start-reverse'},defs);svg('path',{d:'M 0 0 L 10 5 L 0 10 z',fill:'#afbecd'},marker);
 const positions=new Map(),nodes=topology.nodes;const rows=Math.ceil(nodes.length/4);graph.setAttribute('viewBox',`0 0 1160 ${rows*135+50}`);
 nodes.forEach((n,i)=>{const row=Math.floor(i/4),col=row%2?3-i%4:i%4;positions.set(n.id,{x:35+col*287,y:32+row*135});});
 for(const edge of topology.edges){const a=positions.get(edge.source),b=positions.get(edge.target);if(!a||!b)continue;
  let path;if(a.y===b.y){const forward=b.x>a.x,x1=a.x+(forward?230:0),x2=b.x+(forward?0:230);path=`M${x1},${a.y+43} L${x2},${b.y+43}`;}
  else if(a.x===b.x){path=`M${a.x+115},${a.y+86} L${b.x+115},${b.y}`;}
  else{const x1=a.x+115,y1=a.y+86,x2=b.x+115,y2=b.y;const bend=(y1+y2)/2;path=`M${x1},${y1} C${x1},${bend} ${x2},${bend} ${x2},${y2}`;}
  const p=svg('path',{d:path,class:'edge'+(edge.conditional?' conditional':''),'marker-end':'url(#arrow)'},graph);svg('title',{},p).textContent=`${edge.source} → ${edge.target}${edge.label?' · '+edge.label:''}`;
 }
 for(const n of nodes){const pos=positions.get(n.id),state=data?.nodes?.[n.id],boundary=n.id.startsWith('__'),status=state?.status||'WAITING';
  const g=svg('g',{transform:`translate(${pos.x},${pos.y})`,class:`node ${status}${selected===n.id?' selected':''}${boundary?' boundary':''}`,role:'button',tabindex:0,'data-node':n.id,'aria-label':`${n.id}: ${status}`},graph);
  svg('rect',{width:230,height:86,rx:12},g);svg('text',{x:15,y:27,class:'node-name'},g).textContent=n.id;
  svg('text',{x:15,y:49,class:'node-status'},g).textContent=boundary?'Graph boundary':({COMPLETED:'✓ ',RUNNING:'● ',FAILED:'× ',SKIPPED:'− '}[status]||'○ ')+status;
  svg('text',{x:15,y:69,class:'node-time'},g).textContent=boundary?'':`${duration(state?.duration_ms)}${state?.attempt>1?' · attempt '+state.attempt:''}`;
  if(state?.error)svg('title',{},g).textContent=state.error;
  const choose=()=>{selected=n.id;draw();details();};g.onclick=choose;g.onkeydown=e=>{if(e.key==='Enter'||e.key===' '){e.preventDefault();choose();}};
 }
}
function details(){
 const row=data?.nodes?.[selected];$('selectedName').textContent=selected||'노드를 선택하세요';$('detailSummary').replaceChildren();$('detailBody').replaceChildren();
 if(!row){$('detailBody').append(text('p','아직 실행되지 않은 노드입니다.','hint'));return;}
 for(const[label,value]of [['Status',row.status],['Start',time(row.started_at)],['End',time(row.ended_at)],['Duration',duration(row.duration_ms)],['Attempt',row.attempt]]){const cell=text('div','');cell.append(text('small',label),text('span',value));$('detailSummary').append(cell);}
 for(const[key,label]of [['input','Input · bounded preview'],['output','Output · bounded preview'],['state','State · node entry checkpoint'],['metadata','Metadata'],['error','Error']]){const box=document.createElement('details');box.open=key==='output'||key==='error'&&!!row.error;box.append(text('summary',label),text('pre',JSON.stringify(row[key]??null,null,2)));$('detailBody').append(box);}
}
function renderLogs(){const box=$('logs');box.replaceChildren();for(const event of data?.events||[]){const row=text('div','','log-row'+(event.type==='node_error'?' failed':''));row.append(text('time',time(event.timestamp)),text('span',`${event.node||'Run'} · ${event.type}${event.message?'\n'+event.message:''}${event.error?'\n'+event.error:''}`));box.append(row);}if($('follow').checked)box.scrollTop=box.scrollHeight;}
function receive(snapshot){
 const key=`${snapshot.run_id}:${snapshot.version}`;if(lastKey===key)return;lastKey=key;data=snapshot;
 if(snapshot.graph)topology=snapshot.graph;
 if(lastRun!==snapshot.run_id){lastRun=snapshot.run_id;selected=null;history();}
 $('runId').textContent=snapshot.run_id||'아직 실행 기록이 없습니다';$('runStatus').textContent=snapshot.status;
 $('runStatus').style.color=snapshot.status==='FAILED'?'#bd4a4a':snapshot.status==='RUNNING'?'#287be0':'#078477';
 const live=!$('runSelect').value;const running=snapshot.status==='RUNNING';$('start').disabled=requesting||running;$('stop').hidden=!live||!running;
 const active=Object.values(snapshot.nodes||{}).filter(n=>n.status==='RUNNING');$('current').replaceChildren();
 if(active.length){for(const n of active){const row=text('div','','active-item');row.append(text('strong',n.node),text('span',`● RUNNING · ${time(n.started_at)}`));$('current').append(row);}if(!selected)selected=active[0].node;}
 else $('current').append(text('p',snapshot.status==='PAUSED'?'사용자 선택 대기 중':running?'다음 노드 대기 중':snapshot.status==='IDLE'?'실행 대기 중':`실행 종료 · ${snapshot.status}`));
 draw();details();renderLogs();if(snapshot.status!=='RUNNING')history();clock();
}
function clock(){for(const n of document.querySelectorAll('.node.RUNNING')){const row=data?.nodes?.[n.dataset.node];if(row?.started_at)n.querySelector('.node-time').textContent=duration(Date.now()-Date.parse(row.started_at))+(' · attempt '+row.attempt);}
if(!data?.started_at){$('elapsed').textContent='00:00';return;}const end=data.ended_at?Date.parse(data.ended_at):Date.now(),seconds=Math.max(0,Math.floor((end-Date.parse(data.started_at))/1000));$('elapsed').textContent=`${Math.floor(seconds/60).toString().padStart(2,'0')}:${(seconds%60).toString().padStart(2,'0')}`;}
function connect(){
 const epoch=++connectionEpoch;if(transport)transport.close();lastKey='';connection('연결 중…');
 const query=$('runSelect').value?'?run_id='+encodeURIComponent($('runSelect').value):'';
 let received=false;
 const ws=new WebSocket(`${location.protocol==='https:'?'wss:':'ws:'}//${location.host}/ws/workflow${query}`);transport=ws;
 function sse(){if(epoch!==connectionEpoch)return;ws.onclose=null;ws.close();connection('SSE 연결 중…');const source=new EventSource('/api/workflow/events'+query);transport=source;
  source.addEventListener('snapshot',event=>{if(epoch!==connectionEpoch)return;try{receive(JSON.parse(event.data));connection('● Live · SSE');}catch(error){connection('이벤트 읽기 오류',true);}});
  source.onerror=()=>{if(epoch===connectionEpoch)connection('연결 복구 중 · 마지막 상태 유지',true);};
 }
 ws.onmessage=event=>{if(epoch!==connectionEpoch)return;received=true;try{receive(JSON.parse(event.data));connection('● Live · WebSocket');}catch(error){connection('이벤트 읽기 오류',true);}};
 ws.onclose=()=>{if(epoch===connectionEpoch)sse();};ws.onerror=()=>{};
 setTimeout(()=>{if(epoch===connectionEpoch&&!received&&transport===ws)sse();},2500);
}
async function history(){try{const selected=$('runSelect').value,rows=await api('/api/workflow/runs');$('runSelect').replaceChildren(new Option('현재 실행 · Live',''));for(const row of rows)$('runSelect').append(new Option(`${row.run_id} · ${row.status}`,row.run_id));$('runSelect').value=selected;}catch(error){$('runMessage').textContent='이력 로드 실패: '+error.message;}}
$('runSelect').onchange=()=>{selected=null;connect();};
$('runForm').onsubmit=async event=>{event.preventDefault();requesting=true;$('start').disabled=true;$('runMessage').textContent='실행 요청 중…';try{const r=await api('/api/run',{question:$('question').value.trim(),design_mode:'all_eligible',structural_backend:'af3'});$('runMessage').textContent=r.status==='already_running'?'기존 실행이 진행 중입니다.':'파이프라인을 시작했습니다.';$('runSelect').value='';connect();}catch(error){$('runMessage').textContent=error.message;}finally{requesting=false;$('start').disabled=data?.status==='RUNNING';}};
$('stop').onclick=async()=>{$('stop').disabled=true;try{await api('/api/stop',{run_id:data.execution_id});$('runMessage').textContent='중단 요청을 전달했습니다.';}catch(error){$('runMessage').textContent=error.message;}finally{$('stop').disabled=false;}};
(async()=>{try{topology=await api('/api/workflow/graph');draw();await history();connect();}catch(error){connection('연결 실패',true);$('runMessage').textContent=error.message;}})();setInterval(clock,1000);
})();
