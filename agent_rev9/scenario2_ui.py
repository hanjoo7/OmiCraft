"""
시나리오 2 UI — 자원 인지형 후보물질 스크리닝
/scenario2 엔드포인트로 서빙
"""

SCENARIO2_HTML = r"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>OmiCraft — 시나리오 2 · 자원 인지형 후보물질 스크리닝</title>
<link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500&display=swap" rel="stylesheet">
<style>
@import url("https://cdn.jsdelivr.net/gh/orioncactus/pretendard@v1.3.9/dist/web/variable/pretendardvariable-dynamic-subset.min.css");
:root{
  --surface:#F7F9FB;--surface-low:#F2F4F6;--surface-c:#ECEEF0;--surface-high:#E6E8EA;
  --panel:#FFFFFF;--ink:#191C1E;--ink-2:#44474F;--ink-3:#757780;
  --rule:#E0E3E5;--rule-2:#C4C6D0;
  --navy:#001C46;--navy-2:#0B2452;--navy-soft:#D8E2FF;
  --teal:#006A66;--teal-soft:#DFF5F3;--teal-line:#7FCFCB;--teal-txt:#00504D;
  --amber:#9A5B06;--amber-bg:#FDF0D5;--amber-line:#E0B36A;
  --err:#BA1A1A;--err-bg:#FFDAD6;--err-txt:#93000A;
  --dead:#757780;--dead-bg:#E6E8EA;
  --sans:"Pretendard Variable",Inter,system-ui,sans-serif;
  --mono:"JetBrains Mono",monospace;
  --side:260px;
}
*{box-sizing:border-box;margin:0;padding:0}
body{font-family:var(--sans);font-size:15px;line-height:1.6;color:var(--ink);background:var(--surface)}
.app{display:flex;min-height:100vh}
.side{width:var(--side);flex:0 0 var(--side);position:fixed;top:0;bottom:0;left:0;background:var(--panel);border-right:1px solid var(--rule);display:flex;flex-direction:column;padding:20px 0 16px;overflow-y:auto;z-index:20}
.brand{padding:0 20px 16px;display:flex;align-items:center;gap:10px}
.brand b{font-size:19px;font-weight:700;letter-spacing:-.02em}
.cta{margin:0 16px 16px;padding:11px;border:0;border-radius:6px;background:var(--navy);color:#fff;font-family:var(--mono);font-size:13px;font-weight:500;cursor:pointer;text-align:center}
.cta:hover{background:var(--navy-2)}.cta:disabled{opacity:.4;cursor:not-allowed}
.nav{display:flex;flex-direction:column;gap:1px}
.nv{display:flex;align-items:center;gap:10px;padding:10px 20px;font-size:13px;color:var(--ink-2);border-left:3px solid transparent}
.nv.on{background:var(--surface-low);color:var(--teal);font-weight:600;border-left-color:var(--teal)}
.nv.done{color:var(--teal)}.nv .st{margin-left:auto;font-family:var(--mono);font-size:9px;letter-spacing:.06em;padding:2px 7px;border-radius:8px;background:var(--surface-c);color:var(--ink-3)}
.nv.on .st,.nv.done .st{background:var(--teal-soft);color:var(--teal)}
.sres{margin:auto 16px 12px;border:1px solid var(--rule);border-radius:8px;padding:12px 14px}
.srh{font-family:var(--mono);font-size:9px;letter-spacing:.1em;color:var(--ink-3);margin-bottom:10px}
.srb{margin-bottom:8px}.srb .l{display:flex;justify-content:space-between;font-size:11px;margin-bottom:3px}
.srb .l b{font-family:var(--mono);font-size:10px;color:var(--ink)}
.srb .t{height:4px;background:var(--surface-high);border-radius:2px;overflow:hidden}
.srb .t i{display:block;height:100%;width:0;background:var(--teal);border-radius:2px;transition:width .5s}
.main{flex:1;min-width:0;margin-left:var(--side)}
.topbar{position:sticky;top:0;z-index:15;background:rgba(247,249,251,.9);backdrop-filter:blur(12px);border-bottom:1px solid var(--rule);padding:14px 28px;display:flex;align-items:center;gap:16px}
.topbar h1{font-size:20px}.topbar .sub{font-family:var(--mono);font-size:10px;letter-spacing:.07em;color:var(--ink-3);margin-top:2px}
.chip{display:inline-flex;align-items:center;gap:6px;padding:5px 12px;border-radius:999px;border:1px solid var(--rule-2);background:var(--panel);font-family:var(--mono);font-size:11px;color:var(--ink-2)}
.chip i{width:7px;height:7px;border-radius:50%;background:var(--ink-3);flex:0 0 7px}
.chip.run{border-color:var(--teal-line);background:var(--teal-soft);color:var(--teal-txt)}.chip.run i{background:var(--teal);animation:pulse 1.4s infinite}
.chip.ok{border-color:var(--teal-line);background:var(--teal-soft);color:var(--teal-txt)}.chip.ok i{background:var(--teal)}
@keyframes pulse{0%,100%{opacity:1}50%{opacity:.2}}
.content{padding:24px 28px 40px;max-width:1400px}
section{margin-bottom:28px}
.shead{display:flex;align-items:baseline;gap:12px;margin-bottom:12px}
.shead .num{font-family:var(--mono);font-size:11px;color:var(--teal);letter-spacing:.09em;font-weight:500}
.shead h2{font-size:16px}.shead .hint{margin-left:auto;font-size:12px;color:var(--ink-3)}
.card{background:var(--panel);border:1px solid var(--rule);border-radius:8px}
.ch{display:flex;align-items:center;gap:8px;padding:10px 14px;background:var(--surface-low);border-bottom:1px solid var(--rule);border-radius:8px 8px 0 0;font-family:var(--mono);font-size:10px;font-weight:500;letter-spacing:.08em;color:var(--ink-2)}
.ch .r{margin-left:auto;font-size:10px;color:var(--ink-3)}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px;margin-bottom:14px}
.tile{background:var(--panel);border:1px solid var(--rule);border-radius:8px;padding:14px 16px}
.tile .tl{font-size:12px;color:var(--ink-2);margin-bottom:2px}.tile .tv{font-family:var(--mono);font-size:24px;font-weight:500}
/* profile cards */
.profiles{display:grid;grid-template-columns:repeat(3,1fr);gap:12px}
.pcard{background:var(--panel);border:1px solid var(--rule);border-radius:8px;padding:16px;position:relative}
.pcard.best{border-color:var(--teal);box-shadow:0 0 0 2px var(--teal-soft)}
.pcard .pn{font-size:15px;font-weight:600;margin-bottom:4px}.pcard .pt{font-size:12px;color:var(--ink-3);margin-bottom:10px}
.pcard dl{display:grid;grid-template-columns:auto 1fr;gap:4px 12px;font-size:12px;margin:0}
.pcard dt{color:var(--ink-3)}.pcard dd{margin:0;font-family:var(--mono);font-size:11px}
.pcard .badge{position:absolute;top:12px;right:12px;font-family:var(--mono);font-size:9px;padding:3px 8px;border-radius:8px}
.pcard .badge.rec{background:var(--teal-soft);color:var(--teal)}
.pcard .badge.warn{background:var(--amber-bg);color:var(--amber)}
.pcard .badge.danger{background:var(--err-bg);color:var(--err-txt)}
/* cascade table */
.ctbl{width:100%;border-collapse:collapse;font-size:12.5px}
.ctbl th{text-align:left;padding:7px 10px;background:var(--surface-low);border-bottom:1px solid var(--rule);font-family:var(--mono);font-size:10px;color:var(--ink-3)}
.ctbl td{padding:7px 10px;border-bottom:1px solid var(--rule)}
.ctbl .on{color:var(--teal);font-weight:600}.ctbl .skip{color:var(--dead);text-decoration:line-through}
/* log */
.logbox{max-height:350px;overflow-y:auto;padding:12px 16px;font-family:var(--mono);font-size:11.5px;line-height:1.7}
.logbox .le{margin-bottom:2px}.logbox .ts{color:var(--ink-3);margin-right:6px}
.logbox .screening{color:var(--teal)}.logbox .critic{color:var(--err)}.logbox .dossier{color:var(--amber)}
.locked{padding:40px;text-align:center;color:var(--ink-3);font-size:14px}
.foot{padding:20px 28px;margin-left:var(--side);font-size:11px;color:var(--ink-3);line-height:1.7;border-top:1px solid var(--rule)}
</style>
</head>
<body>
<div class="app">
<div class="side">
  <div class="brand"><b>OmiCraft</b></div>
  <button class="cta" id="btnRun" onclick="go()">▶ 스크리닝 실행</button>
  <div class="nav" id="agentNav"></div>
  <div class="sres">
    <div class="srh">RESOURCE USAGE</div>
    <div class="srb"><div class="l"><span>Elapsed</span><b id="rvT">0s</b></div><div class="t"><i id="riT"></i></div></div>
    <div class="srb"><div class="l"><span>Evidence Cards</span><b id="rvE">0</b></div><div class="t"><i id="riE"></i></div></div>
    <div class="srb"><div class="l"><span>Skills Used</span><b id="rvS">0/4</b></div><div class="t"><i id="riS"></i></div></div>
  </div>
</div>
<div class="main">
<div class="topbar">
  <div><h1>자원 인지형 후보물질 스크리닝</h1><div class="sub">시나리오 2</div></div>
  <span style="flex:1"></span>
  <span class="chip" id="statusChip"><i></i><span id="statusText">IDLE</span></span>
</div>
<div class="content">

<section>
  <div class="shead"><span class="num">01</span><h2>입력과 자원 프로파일</h2><span class="hint">같은 표적·같은 라이브러리 — 자원만 바꾸면 시스템이 다른 경로를 선택합니다</span></div>
  <div class="tiles">
    <div class="tile"><div class="tl">입력 표적</div><div class="tv" id="tvTarget">—</div></div>
    <div class="tile"><div class="tl">라이브러리</div><div class="tv">240K</div></div>
    <div class="tile"><div class="tl">Skills</div><div class="tv" id="tvSkills">—</div></div>
  </div>
  <div class="profiles" id="profileCards"><div class="locked">실행 후 표시됩니다.</div></div>
</section>

<section>
  <div class="shead"><span class="num">02</span><h2>자원 인지 캐스케이드</h2><span class="hint">프로파일별 단계 ON/SKIP</span></div>
  <div class="card"><div id="cascadeBox" style="padding:16px"><div class="locked">실행 후 표시됩니다.</div></div></div>
</section>

<section>
  <div class="shead"><span class="num">03</span><h2>프로파일 대조</h2></div>
  <div class="card"><div id="compareBox" style="padding:16px"><div class="locked">실행 후 표시됩니다.</div></div></div>
</section>

<section>
  <div class="shead"><span class="num">04</span><h2>불일치 탐지</h2><span class="hint">docking · Boltz-2 결과가 충돌하면 Critic이 개입합니다</span></div>
  <div class="card"><div id="disagreeBox" style="padding:16px"><div class="locked">실행 후 표시됩니다.</div></div></div>
</section>

<section>
  <div class="shead"><span class="num">05</span><h2>실행 로그</h2><span class="hint" id="logCount">0건</span></div>
  <div class="card"><div class="logbox" id="logBox"><div class="locked">실행 후 표시됩니다.</div></div></div>
</section>

</div>
<div class="foot">OmiCraft · 시나리오 2 · OMICS2DRUG · 실제 에이전트 실행 결과</div>
</div>
</div>

<script>
const $=id=>document.getElementById(id);
const AGENTS=[
  {id:'planner',n:'Planner'},{id:'discovery',n:'Discovery'},{id:'qualification',n:'Qualification'},
  {id:'design',n:'Design'},{id:'critic',n:'Critic'},{id:'screening',n:'Screening'},{id:'dossier_s2',n:'Dossier'}
];
let ws,ti,logN=0;

$('agentNav').innerHTML=AGENTS.map(a=>`<div class="nv" id="nav-${a.id}"><span>${a.n}</span><span class="st" id="st-${a.id}">IDLE</span></div>`).join('');

function addLog(node,msg){
  if(logN===0)$('logBox').innerHTML='';
  const cls=node==='screening'?'screening':node==='critic'?'critic':node.includes('dossier')?'dossier':'';
  const div=document.createElement('div');div.className='le';
  div.innerHTML=`<span class="ts ${cls}">[${node}]</span> ${esc(msg)}`;
  $('logBox').appendChild(div);$('logBox').scrollTop=$('logBox').scrollHeight;
  logN++;$('logCount').textContent=logN+'건';
}
function esc(s){const d=document.createElement('span');d.textContent=s;return d.innerHTML;}

function buildProfiles(target){
  const profiles=[
    {id:'A',name:'80GB GPU × 1',vram:'80 GB',gpu:'6.0 h',batch:'대형 가능',badge:'rec',badgeText:'RECOMMENDED'},
    {id:'B',name:'48GB GPU × 2',vram:'48 GB × 2',gpu:'6.0 h',batch:'축소 필요',badge:'warn',badgeText:'REDUCED'},
    {id:'C',name:'GPU 없음 (CPU only)',vram:'0',gpu:'0 h',batch:'CPU only',badge:'danger',badgeText:'TIER 3'},
  ];
  $('profileCards').innerHTML=profiles.map(p=>`
    <div class="pcard${p.id==='A'?' best':''}">
      <div class="pn">${p.name}</div>
      <div class="pt">프로파일 ${p.id}</div>
      <dl><dt>VRAM</dt><dd>${p.vram}</dd><dt>GPU 예산</dt><dd>${p.gpu}</dd><dt>배치</dt><dd>${p.batch}</dd></dl>
      <span class="badge ${p.badge}">${p.badgeText}</span>
    </div>
  `).join('');
}

function buildCascade(screening){
  if(!screening||!screening.profiles)return;
  const pids=['A','B','C'];
  let html='<table class="ctbl"><thead><tr><th>단계</th><th>도구</th>';
  pids.forEach(p=>html+=`<th>P${p}</th>`);
  html+='</tr></thead><tbody>';
  const refCascade=screening.profiles.A?.cascade||[];
  refCascade.forEach((step,i)=>{
    html+='<tr>';
    html+=`<td>${step.name}</td><td style="font-size:11px">${step.tool}</td>`;
    pids.forEach(pid=>{
      const ps=screening.profiles[pid]?.cascade?.[i];
      if(ps&&ps.status==='ON'){
        html+=`<td class="on">${ps.input.toLocaleString()} → ${ps.pass.toLocaleString()}<br><span style="font-size:10px;color:var(--ink-3)">GPU ${ps.gpu_h}h</span></td>`;
      } else {
        html+=`<td class="skip">SKIP</td>`;
      }
    });
    html+='</tr>';
  });
  html+='</tbody></table>';
  $('cascadeBox').innerHTML=html;
}

function buildCompare(screening){
  if(!screening||!screening.profiles)return;
  let html='<table class="ctbl"><thead><tr><th>프로파일</th><th>Virtual Hits</th><th>GPU</th><th>단계</th><th>불일치</th><th>Tier</th></tr></thead><tbody>';
  ['A','B','C'].forEach(pid=>{
    const r=screening.profiles[pid];if(!r)return;
    const tier=r.all_tier3?'<span style="color:var(--err)">Tier 3</span>':'<span style="color:var(--teal)">Tier 1-2</span>';
    const best=pid==='A'?' style="background:var(--teal-soft)"':'';
    html+=`<tr${best}><td style="font-weight:600">프로파일 ${pid}</td>`;
    html+=`<td style="font-family:var(--mono);font-size:16px;font-weight:600">${r.final_virtual_hits}</td>`;
    html+=`<td>${r.total_gpu_h}h</td>`;
    html+=`<td>${r.n_active_steps}/${r.n_active_steps+r.n_skipped_steps}</td>`;
    html+=`<td>${r.disagreements?.length||0}건</td>`;
    html+=`<td>${tier}</td></tr>`;
  });
  html+='</tbody></table>';
  $('compareBox').innerHTML=html;
}

function buildDisagreements(screening){
  if(!screening||!screening.profiles)return;
  const allD=[];
  ['A','B','C'].forEach(pid=>{
    (screening.profiles[pid]?.disagreements||[]).forEach(d=>{
      allD.push({...d,profile:pid});
    });
  });
  if(!allD.length){$('disagreeBox').innerHTML='<div style="padding:20px;color:var(--ink-3)">불일치 없음</div>';return;}
  let html='<table class="ctbl"><thead><tr><th>프로파일</th><th>화합물</th><th>이슈</th><th>해소</th></tr></thead><tbody>';
  allD.forEach(d=>{
    const resolved=d.resolved?'<span style="color:var(--teal)">✓ 해소</span>':'<span style="color:var(--err)">✗ 미해결</span>';
    html+=`<tr><td>P${d.profile}</td><td style="font-family:var(--mono)">${d.compound}</td><td style="font-size:11px">${d.issue}</td><td>${resolved}<br><span style="font-size:10px;color:var(--ink-3)">${d.action}</span></td></tr>`;
  });
  html+='</tbody></table>';
  $('disagreeBox').innerHTML=html;
}

function conn(){
  ws=new WebSocket('ws://'+location.host+'/ws');
  ws.onmessage=e=>handle(JSON.parse(e.data));
  ws.onclose=()=>setTimeout(conn,2000);
}

function handle(d){
  if(d.type==='node'){
    AGENTS.forEach(a=>{const el=$('nav-'+a.id);if(el&&el.classList.contains('on')){el.className='nv done';$('st-'+a.id).textContent='DONE';}});
    const nav=$('nav-'+d.node);if(nav){nav.className='nv on';$('st-'+d.node).textContent='RUNNING';}
    $('rvT').textContent=d.elapsed+'s';$('riT').style.width=Math.min(100,d.elapsed/180*100)+'%';
    (d.messages||[]).forEach(m=>addLog(d.node,m));

    // Extract screening data from messages
    (d.messages||[]).forEach(m=>{
      let match;
      if(match=m.match(/스크리닝 대상 표적:\s*(\S+)/)) $('tvTarget').textContent=match[1];
      if(match=m.match(/Skill.*?:\s*[✓✗]/)) {
        const n=document.querySelectorAll('.le').length;
        const skills=(m.match(/✓/g)||[]).length;
        $('tvSkills').textContent=skills+'/4';$('rvS').textContent=skills+'/4';$('riS').style.width=skills/4*100+'%';
      }
    });

    if(d.verdict){
      addLog('critic','판정: '+d.verdict);
      $('statusChip').className='chip '+(d.verdict==='APPROVE'?'ok':'run');
      $('statusText').textContent=d.verdict;
    }
    $('logBox').scrollTop=$('logBox').scrollHeight;
  }
  if(d.type==='done'){
    clearInterval(ti);
    AGENTS.forEach(a=>{$('nav-'+a.id).className='nv done';$('st-'+a.id).textContent='DONE';});
    $('statusChip').className='chip ok';$('statusText').textContent=d.final_verdict||'DONE';
    $('rvT').textContent=d.total_time+'s';$('btnRun').disabled=false;

    // Load screening results
    fetch('/api/screening').then(r=>r.json()).then(data=>{
      if(data.target)$('tvTarget').textContent=data.target;
      buildProfiles(data.target);
      buildCascade(data);
      buildCompare(data);
      buildDisagreements(data);
    }).catch(()=>{});
  }
}

async function go(){
  $('btnRun').disabled=true;$('statusChip').className='chip run';$('statusText').textContent='RUNNING';
  logN=0;$('logBox').innerHTML='';
  ['cascadeBox','compareBox','disagreeBox'].forEach(id=>$(id).innerHTML='<div class="locked">실행 중...</div>');
  $('profileCards').innerHTML='<div class="locked">실행 중...</div>';
  ['tvTarget','tvSkills'].forEach(id=>$(id).textContent='—');
  AGENTS.forEach(a=>{$('nav-'+a.id).className='nv';$('st-'+a.id).textContent='IDLE';});
  $('rvT').textContent='0s';$('rvE').textContent='0';$('rvS').textContent='0/4';
  ['riT','riE','riS'].forEach(id=>$(id).style.width='0');

  const st=Date.now();
  ti=setInterval(()=>{$('rvT').textContent=((Date.now()-st)/1000).toFixed(1)+'s';$('riT').style.width=Math.min(100,(Date.now()-st)/180000*100)+'%';},200);

  await fetch('/api/run_s2',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({question:'TNBC 표적 발굴 + 자원 인지형 스크리닝'})});
}

conn();
</script>
</body>
</html>"""
