"""
시나리오 1 스타일 UI — 실제 에이전트 실행 결과를 시나리오 1 데모와 동일한 UI로 표시
/scenario1 엔드포인트로 서빙
"""

SCENARIO1_HTML = r"""<!DOCTYPE html>
<html lang="ko">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>OmiCraft — 시나리오 1 · 표적 발굴에서 모달리티 선택까지</title>
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
.brand span{font-family:var(--mono);font-size:10px;color:var(--ink-3);letter-spacing:.06em}
.cta{margin:0 16px 16px;padding:11px;border:0;border-radius:6px;background:var(--navy);color:#fff;font-family:var(--mono);font-size:13px;font-weight:500;cursor:pointer;text-align:center}
.cta:hover{background:var(--navy-2)}.cta:disabled{opacity:.4;cursor:not-allowed}

.nav{display:flex;flex-direction:column;gap:1px}
.nv{display:flex;align-items:center;gap:10px;padding:10px 20px;font-size:13px;color:var(--ink-2);border-left:3px solid transparent;cursor:default}
.nv.on{background:var(--surface-low);color:var(--teal);font-weight:600;border-left-color:var(--teal)}
.nv .st{margin-left:auto;font-family:var(--mono);font-size:9px;letter-spacing:.06em;padding:2px 7px;border-radius:8px;background:var(--surface-c);color:var(--ink-3)}
.nv.on .st{background:var(--teal-soft);color:var(--teal)}
.nv.done .st{background:var(--teal-soft);color:var(--teal)}

.sres{margin:auto 16px 12px;border:1px solid var(--rule);border-radius:8px;padding:12px 14px}
.srh{font-family:var(--mono);font-size:9px;letter-spacing:.1em;color:var(--ink-3);margin-bottom:10px}
.srb{margin-bottom:8px}.srb .l{display:flex;justify-content:space-between;font-size:11px;color:var(--ink-2);margin-bottom:3px}
.srb .l b{font-family:var(--mono);font-size:10px;color:var(--ink)}
.srb .t{height:4px;background:var(--surface-high);border-radius:2px;overflow:hidden}
.srb .t i{display:block;height:100%;width:0;background:var(--teal);border-radius:2px;transition:width .5s}

.main{flex:1;min-width:0;margin-left:var(--side)}
.topbar{position:sticky;top:0;z-index:15;background:rgba(247,249,251,.9);backdrop-filter:blur(12px);border-bottom:1px solid var(--rule);padding:14px 28px;display:flex;align-items:center;gap:16px}
.topbar h1{font-size:20px}
.topbar .sub{font-family:var(--mono);font-size:10px;letter-spacing:.07em;color:var(--ink-3);margin-top:2px}
.chip{display:inline-flex;align-items:center;gap:6px;padding:5px 12px;border-radius:999px;border:1px solid var(--rule-2);background:var(--panel);font-family:var(--mono);font-size:11px;color:var(--ink-2)}
.chip i{width:7px;height:7px;border-radius:50%;background:var(--ink-3);flex:0 0 7px}
.chip.run{border-color:var(--teal-line);background:var(--teal-soft);color:var(--teal-txt)}
.chip.run i{background:var(--teal);animation:pulse 1.4s infinite}
.chip.ok{border-color:var(--teal-line);background:var(--teal-soft);color:var(--teal-txt)}
.chip.ok i{background:var(--teal)}
@keyframes pulse{0%,100%{opacity:1}50%{opacity:.2}}

.content{padding:24px 28px 40px;max-width:1400px}
section{margin-bottom:28px}
.shead{display:flex;align-items:baseline;gap:12px;margin-bottom:12px}
.shead .num{font-family:var(--mono);font-size:11px;color:var(--teal);letter-spacing:.09em;font-weight:500}
.shead h2{font-size:16px}
.shead .hint{margin-left:auto;font-size:12px;color:var(--ink-3)}
.card{background:var(--panel);border:1px solid var(--rule);border-radius:8px}
.ch{display:flex;align-items:center;gap:8px;padding:10px 14px;background:var(--surface-low);border-bottom:1px solid var(--rule);border-radius:8px 8px 0 0;font-family:var(--mono);font-size:10px;font-weight:500;letter-spacing:.08em;color:var(--ink-2)}
.ch .r{margin-left:auto;font-size:10px;color:var(--ink-3)}

.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(160px,1fr));gap:10px;margin-bottom:14px}
.tile{background:var(--panel);border:1px solid var(--rule);border-radius:8px;padding:14px 16px}
.tile .tl{font-size:12px;color:var(--ink-2);margin-bottom:2px}
.tile .tv{font-family:var(--mono);font-size:24px;font-weight:500;letter-spacing:-.02em}

/* stepper */
.stepper{display:flex;padding:16px 20px 14px;gap:0;align-items:flex-start}
.sp{display:flex;flex-direction:column;align-items:center;gap:7px;flex:1;text-align:center}
.sp .dot{width:28px;height:28px;border-radius:50%;border:1.5px solid var(--rule-2);background:var(--panel);display:grid;place-items:center;font-family:var(--mono);font-size:10px;color:var(--ink-3);transition:.3s}
.sp .lb{font-size:11px;color:var(--ink-3);transition:.3s}
.sp.act .dot{border-color:var(--teal);color:var(--teal);box-shadow:0 0 0 3px var(--teal-soft)}
.sp.act .lb{color:var(--teal);font-weight:600}
.sp.done .dot{background:var(--teal);border-color:var(--teal);color:#fff}
.sp.done .lb{color:var(--teal-txt)}
.sp .line{width:100%;height:2px;background:var(--rule);margin-top:13px}
.sp.done .line{background:var(--teal)}

/* funnel */
.funnel{padding:18px 20px}
.fb{display:flex;align-items:center;gap:12px;margin-bottom:6px}
.fb .bar{height:22px;border-radius:4px;background:var(--navy-soft);transition:width .6s;min-width:4px;display:flex;align-items:center;padding-left:8px;font-family:var(--mono);font-size:10px;color:var(--navy);font-weight:500}
.fb.cut .bar{background:var(--err-bg);color:var(--err-txt)}
.fb.final .bar{background:var(--teal-soft);color:var(--teal-txt)}
.fb .fl{font-size:11.5px;color:var(--ink-2);min-width:200px;text-align:right}

/* candidate table */
.ctbl{width:100%;border-collapse:collapse;font-size:13px}
.ctbl th{text-align:left;padding:8px 12px;background:var(--surface-low);border-bottom:1px solid var(--rule);font-family:var(--mono);font-size:10px;letter-spacing:.06em;color:var(--ink-3);font-weight:500}
.ctbl td{padding:8px 12px;border-bottom:1px solid var(--rule)}
.ctbl .gene{font-weight:600;color:var(--navy)}
.vbadge{display:inline-block;font-weight:600;padding:2px 8px;border-radius:4px;font-size:11px;font-family:var(--mono)}
.vbadge.adv{background:var(--teal-soft);color:var(--teal)}
.vbadge.rej{background:var(--err-bg);color:var(--err-txt)}
.vbadge.hold{background:var(--amber-bg);color:var(--amber)}

/* modality routing */
.rtbl{width:100%;border-collapse:collapse;font-size:12.5px}
.rtbl th{text-align:left;padding:7px 12px;background:var(--surface-low);border-bottom:1px solid var(--rule);font-size:10px;font-family:var(--mono);color:var(--ink-3)}
.rtbl td{padding:7px 12px;border-bottom:1px solid var(--rule)}
.rtbl .hl{background:var(--teal-soft)}

/* dossier */
.dos{border:1px solid var(--teal-line);border-radius:8px;margin-bottom:14px;overflow:hidden}
.dos .dh{padding:12px 16px;background:var(--teal-soft);display:flex;align-items:center;gap:10px}
.dos .dh b{font-size:15px;color:var(--teal-txt)}.dos .dh .tag{margin-left:auto;font-family:var(--mono);font-size:10px;color:var(--teal);letter-spacing:.06em}
.dos .db{padding:14px 16px}
.dos dl{display:grid;grid-template-columns:140px 1fr;gap:6px 14px;font-size:12.5px;margin:0}
.dos dt{font-weight:600;color:var(--ink-2)}.dos dd{margin:0}
.dos .wet{margin-top:10px;padding-top:10px;border-top:1px solid var(--rule)}
.dos .wet li{font-size:12px;margin-bottom:3px;color:var(--ink-2)}

/* log */
.logbox{max-height:300px;overflow-y:auto;padding:12px 16px;font-family:var(--mono);font-size:11.5px;line-height:1.7}
.logbox .le{margin-bottom:2px}
.logbox .ts{color:var(--ink-3);margin-right:6px}
.logbox .plan{color:#8957e5}.logbox .omics{color:#2f81f7}.logbox .qual{color:var(--teal)}.logbox .design{color:var(--amber)}.logbox .crit{color:var(--err)}.logbox .safe{color:var(--err)}

.locked{padding:40px;text-align:center;color:var(--ink-3);font-size:14px}
.foot{padding:20px 28px;margin-left:var(--side);font-size:11px;color:var(--ink-3);line-height:1.7;border-top:1px solid var(--rule)}
</style>
</head>
<body>
<div class="app">

<!-- SIDEBAR -->
<div class="side">
  <div class="brand"><b>OmiCraft</b></div>
  <button class="cta" id="btnRun" onclick="startRun()">▶ 분석 실행</button>

  <div class="nav" id="agentNav"></div>

  <div class="sres">
    <div class="srh">RESOURCE USAGE</div>
    <div class="srb"><div class="l"><span>LLM Calls</span><b id="rvL">0</b></div><div class="t"><i id="riL"></i></div></div>
    <div class="srb"><div class="l"><span>Evidence Cards</span><b id="rvE">0</b></div><div class="t"><i id="riE"></i></div></div>
    <div class="srb"><div class="l"><span>Elapsed</span><b id="rvT">0s</b></div><div class="t"><i id="riT"></i></div></div>
  </div>
</div>

<!-- MAIN -->
<div class="main">
<div class="topbar">
  <div><h1>표적 발굴에서 모달리티 선택까지</h1><div class="sub">시나리오 1 · TNBC</div></div>
  <span style="flex:1"></span>
  <span class="chip" id="statusChip"><i></i><span id="statusText">IDLE</span></span>
</div>

<div class="content">

<!-- 01 연구질문 -->
<section>
  <div class="shead"><span class="num">01</span><h2>연구질문</h2></div>
  <div class="card">
    <div style="padding:18px 20px">
      <div style="font-size:17px;font-weight:600;line-height:1.4" id="questionText">삼중음성 유방암(TNBC)에서 개발 가능한 치료 표적과 그에 적합한 모달리티를 제안하라</div>
      <div style="margin-top:10px;font-size:13px;color:var(--ink-2)" id="planSummary">분석을 실행하면 코호트·비교군 정보가 표시됩니다.</div>
    </div>
  </div>
</section>

<!-- 02 진행 Stepper -->
<section>
  <div class="shead"><span class="num">02</span><h2>실행 진행</h2></div>
  <div class="card"><div class="stepper" id="stepper"></div></div>
</section>

<!-- 03 핵심 지표 -->
<section>
  <div class="shead"><span class="num">03</span><h2>핵심 지표</h2></div>
  <div class="tiles">
    <div class="tile"><div class="tl">DEG (유의 유전자)</div><div class="tv" id="tvDEG">—</div></div>
    <div class="tile"><div class="tl">Cell-of-origin 확인</div><div class="tv" id="tvCO">—</div></div>
    <div class="tile"><div class="tl">ADVANCE</div><div class="tv" id="tvADV">—</div></div>
    <div class="tile"><div class="tl">REJECT (Safety)</div><div class="tv" id="tvREJ">—</div></div>
    <div class="tile"><div class="tl">Evidence Cards</div><div class="tv" id="tvEV">—</div></div>
    <div class="tile"><div class="tl">Critic 판정</div><div class="tv" id="tvCR">—</div></div>
  </div>
</section>

<!-- 03b 오믹스 분석 결과 -->
<section>
  <div class="shead"><span class="num">03b</span><h2>오믹스 분석 결과</h2><span class="hint">Omics Discovery Agent 원 산출물</span></div>
  <div style="display:grid;grid-template-columns:1fr 1fr;gap:12px" id="chartsGrid">

    <div class="card">
      <div class="ch">SAMPLE DISTANCE MATRIX<span class="r" id="rDist">대기</span></div>
      <div style="padding:12px"><canvas id="cvDist" width="480" height="300" style="width:100%;height:auto"></canvas></div>
      <div style="padding:0 14px 10px;font-size:11px;color:var(--ink-3)">비교군 구성 — TNBC vs non-TNBC 샘플 간 거리 히트맵</div>
    </div>

    <div class="card">
      <div class="ch">PCA — SAMPLE CLUSTERING<span class="r" id="rPca">대기</span></div>
      <div style="padding:12px"><canvas id="cvPca" width="480" height="300" style="width:100%;height:auto"></canvas></div>
      <div style="padding:0 14px 10px;font-size:11px;color:var(--ink-3)">주성분 분석 — TNBC(적) vs non-TNBC(청) 군집 구조</div>
    </div>

    <div class="card">
      <div class="ch">VOLCANO — DIFFERENTIAL EXPRESSION<span class="r" id="rVol">대기</span></div>
      <div style="padding:12px"><canvas id="cvVol" width="480" height="300" style="width:100%;height:auto"></canvas></div>
      <div style="padding:0 14px 10px;font-size:11px;color:var(--ink-3)">DESeq2 결과 — x: log2FC, y: -log10(FDR), 빨강: 상향, 청록: 하향</div>
    </div>

    <div class="card">
      <div class="ch">GSEA — ENRICHMENT<span class="r" id="rGsea">대기</span></div>
      <div style="padding:12px"><canvas id="cvGsea" width="480" height="300" style="width:100%;height:auto"></canvas></div>
      <div style="padding:0 14px 10px;font-size:11px;color:var(--ink-3)">MSigDB Hallmark · KEGG · Reactome — NES 상위 유의 경로</div>
    </div>

  </div>
</section>

<!-- 04 후보 Funnel -->
<section>
  <div class="shead"><span class="num">04</span><h2>후보 Funnel</h2><span class="hint">단계별 후보 수 감소</span></div>
  <div class="card"><div class="funnel" id="funnelBox"><div class="locked">분석 실행 후 표시됩니다.</div></div></div>
</section>

<!-- 05 후보 판정 -->
<section>
  <div class="shead"><span class="num">05</span><h2>후보 판정</h2><span class="hint">ADVANCE · REJECT · HOLD</span></div>
  <div class="card"><div id="candBox" style="padding:16px"><div class="locked">분석 실행 후 표시됩니다.</div></div></div>
</section>

<!-- 06 모달리티 라우팅 -->
<section>
  <div class="shead"><span class="num">06</span><h2>모달리티 라우팅</h2></div>
  <div class="card"><div id="routeBox" style="padding:16px"><div class="locked">분석 실행 후 표시됩니다.</div></div></div>
</section>

<!-- 07 실행 로그 -->
<section>
  <div class="shead"><span class="num">07</span><h2>실행 로그</h2><span class="hint" id="logCount">0건</span></div>
  <div class="card"><div class="logbox" id="logBox"><div class="locked">분석 실행 후 표시됩니다.</div></div></div>
</section>

<!-- 08 Dossier -->
<section>
  <div class="shead"><span class="num">08</span><h2>표적 Dossier</h2><span class="hint">근거·등급·반대 근거·검증계획</span></div>
  <div id="dosBox"><div class="locked">분석 완료 후 Dossier가 생성됩니다.</div></div>
</section>

</div>

<div class="foot">
  OmiCraft · 시나리오 1 · OMICS2DRUG · 실제 에이전트 실행 결과를 표시합니다.
</div>
</div>
</div>

<script>
const $=id=>document.getElementById(id);
const AGENTS=[
  {id:'planner',n:'Coordinator–Planner',d:'연구질문 → 실행 DAG'},
  {id:'discovery',n:'Omics Discovery',d:'환자 데이터 → 후보 풀'},
  {id:'qualification',n:'Target Qualification',d:'세포 맥락 → 개발 경로'},
  {id:'design',n:'Therapeutic Design',d:'modality → 설계 명세'},
  {id:'critic',n:'Critic Agent',d:'독립 검증 · 피드백'},
  {id:'dossier',n:'Dossier',d:'최종 산출물'}
];
const STEPS=['Planning','Discovery','Qualification','Design','Critic','Dossier'];
const STEP_MAP={planner:0,discovery:1,qualification:2,design:3,critic:4,dossier:5};

let ws,ti,evN=0,logN=0,advTargets=[];

// build sidebar nav
$('agentNav').innerHTML=AGENTS.map(a=>`<div class="nv" id="nav-${a.id}"><span>${a.n}</span><span class="st" id="st-${a.id}">IDLE</span></div>`).join('');

// build stepper
$('stepper').innerHTML=STEPS.map((s,i)=>`<div class="sp" id="sp${i}"><div class="dot">${i+1}</div><div class="lb">${s}</div></div>`).join('');

function setStep(idx,state){
  for(let k=0;k<STEPS.length;k++){
    const el=$('sp'+k);if(!el)continue;
    el.className='sp'+(k<idx?' done':k===idx?(state==='done'?' done':' act'):'');
  }
}

function addLog(node,msg){
  if(logN===0)$('logBox').innerHTML='';
  const cls=node==='planner'?'plan':node==='discovery'?'omics':node==='qualification'?'qual':node==='design'?'design':node==='critic'?'crit':'';
  const div=document.createElement('div');
  div.className='le';
  div.innerHTML=`<span class="ts ${cls}">[${node}]</span> ${esc(msg)}`;
  $('logBox').appendChild(div);
  $('logBox').scrollTop=$('logBox').scrollHeight;
  logN++;$('logCount').textContent=logN+'건';
}

function esc(s){const d=document.createElement('span');d.textContent=s;return d.innerHTML;}

function buildFunnel(data){
  // data: array of {label, n, cut, final}
  const maxN=Math.max(...data.map(d=>d.n));
  $('funnelBox').innerHTML=data.map(d=>{
    const w=Math.max(4,d.n/maxN*100);
    const cls=d.final?'final':d.cut?'cut':'';
    return `<div class="fb ${cls}"><div class="fl">${d.label}</div><div class="bar" style="width:${w}%">${d.n.toLocaleString()}</div></div>`;
  }).join('');
}

function buildCandidates(targets){
  let html='<table class="ctbl"><thead><tr><th>ID</th><th>GENE</th><th>TIER</th><th>MODALITY</th><th>VERDICT</th></tr></thead><tbody>';
  targets.forEach((t,i)=>{
    const cls=t.judgment==='ADVANCE'?'adv':t.judgment==='REJECT'?'rej':'hold';
    html+=`<tr><td style="font-family:var(--mono);font-size:12px">T-${String(i+1).padStart(2,'0')}</td>`;
    html+=`<td class="gene">${t.gene_name||''}</td>`;
    html+=`<td style="font-family:var(--mono);font-size:12px">${t.bio_tier||''}+${t.dev_tier||''}</td>`;
    html+=`<td>${t.modality||'—'}</td>`;
    html+=`<td><span class="vbadge ${cls}">${t.judgment||''}</span></td></tr>`;
  });
  html+='</tbody></table>';
  $('candBox').innerHTML=html;
}

function buildRouting(targets){
  const routes=[
    {k:'SURFACE',m:'Antibody / ADC',c:'종양 특이 세포 표면 단백질 · internalization'},
    {k:'INTRACELLULAR',m:'Small molecule',c:'세포 내 효소 · pocket 확보'},
    {k:'DEGRADER',m:'Degrader (PROTAC)',c:'세포 내 단백질 · pocket 부적합'},
    {k:'SECRETED',m:'중화항체 / ligand trap',c:'병적 분비 ligand'},
    {k:'RNA_THERAPEUTIC',m:'siRNA / ASO',c:'과발현 RNA'},
    {k:'BIOMARKER',m:'Biomarker · 동반진단',c:'직접 조절 부적합 · 층화 가치'},
  ];
  const advMods={};
  targets.filter(t=>t.judgment==='ADVANCE').forEach(t=>{
    const m=t.modality||'BIOMARKER';
    if(!advMods[m])advMods[m]=[];
    advMods[m].push(t.gene_name);
  });
  let html='<table class="rtbl"><thead><tr><th>후보 유형</th><th>모달리티</th><th>적용 표적</th></tr></thead><tbody>';
  routes.forEach(r=>{
    const genes=advMods[r.k]||advMods[r.m]||[];
    const hl=genes.length>0?' class="hl"':'';
    html+=`<tr${hl}><td>${r.c}</td><td style="font-weight:600">${r.m}</td><td>${genes.join(', ')||'—'}</td></tr>`;
  });
  html+='</tbody></table>';
  $('routeBox').innerHTML=html;
}

function buildDossier(targets){
  const adv=targets.filter(t=>t.judgment==='ADVANCE').slice(0,3);
  if(!adv.length){$('dosBox').innerHTML='<div class="locked">ADVANCE 표적이 없습니다.</div>';return;}
  $('dosBox').innerHTML=adv.map((t,i)=>{
    const wetlab=t.wetlab_plan||['결합 친화도 측정','세포주 패널 검증','정상조직 counter-screen','In vivo 효능 시험'];
    return `<div class="dos">
      <div class="dh"><b>T-${String(i+1).padStart(2,'0')} · ${t.gene_name}</b><span class="tag">ADVANCE · ${t.bio_tier}+${t.dev_tier}</span></div>
      <div class="db">
        <dl>
          <dt>TARGET</dt><dd>${t.gene_name}</dd>
          <dt>BIO CONFIDENCE</dt><dd>${t.bio_tier}</dd>
          <dt>DEV READINESS</dt><dd>${t.dev_tier}</dd>
          <dt>MODALITY</dt><dd>${t.modality}</dd>
          <dt>SAFETY</dt><dd>${t.safety_verdict||'PASS'}</dd>
          <dt>log2FC</dt><dd>${(t.log2fc||0).toFixed(2)}</dd>
        </dl>
        <div class="wet"><b style="font-size:12px">Wet-lab 검증 계획</b><ul>${wetlab.map(w=>'<li>'+w+'</li>').join('')}</ul></div>
      </div>
    </div>`;
  }).join('');
}

// WebSocket
function conn(){
  ws=new WebSocket('ws://'+location.host+'/ws');
  ws.onmessage=e=>handle(JSON.parse(e.data));
  ws.onclose=()=>setTimeout(conn,2000);
}

function handle(d){
  if(d.type==='node'){
    const idx=STEP_MAP[d.node];
    if(idx!==undefined)setStep(idx,'act');

    // sidebar state
    AGENTS.forEach(a=>{const el=$('st-'+a.id);if(el)el.textContent='DONE';$('nav-'+a.id).className='nv done';});
    const navEl=$('nav-'+d.node);if(navEl){navEl.className='nv on';$('st-'+d.node).textContent='RUNNING';}

    // resource
    $('rvT').textContent=d.elapsed+'s';
    $('riT').style.width=Math.min(100,d.elapsed/120*100)+'%';

    // messages → log
    (d.messages||[]).forEach(m=>addLog(d.node,m));

    // verdict
    if(d.verdict){
      addLog('critic','판정: '+d.verdict);
      $('tvCR').textContent=d.verdict;
      $('tvCR').style.color=d.verdict==='APPROVE'?'var(--teal)':d.verdict==='HOLD'?'var(--amber)':'var(--err)';
      $('statusChip').className='chip '+(d.verdict==='APPROVE'?'ok':'run');
      $('statusText').textContent=d.verdict;
    }

    // advance targets
    if(d.advance_targets&&d.advance_targets.length){
      $('tvADV').textContent=d.advance_targets.length;
    }

    // extract KPIs from messages
    (d.messages||[]).forEach(m=>{
      let match;
      if(match=m.match(/DEG:\s*(\d+)\s*significant/)) $('tvDEG').textContent=parseInt(match[1]).toLocaleString();
      if(match=m.match(/TUMOR=(\d+)/)) $('tvCO').textContent=match[1];
      if(match=m.match(/REJECT=(\d+)/)) $('tvREJ').textContent=match[1];
      if(match=m.match(/Evidence Cards?:\s*(\d+)/i)) {evN=parseInt(match[1]);$('tvEV').textContent=evN;$('rvE').textContent=evN;$('riE').style.width=Math.min(100,evN/20*100)+'%';}
      if(match=m.match(/ADVANCE=(\d+)/)) $('tvADV').textContent=match[1];
    });
  }

  if(d.type==='done'){
    clearInterval(ti);
    setStep(STEPS.length-1,'done');
    AGENTS.forEach(a=>{$('nav-'+a.id).className='nv done';$('st-'+a.id).textContent='DONE';});
    $('statusChip').className='chip ok';
    $('statusText').textContent=d.final_verdict||'DONE';
    $('rvT').textContent=d.total_time+'s';
    $('btnRun').disabled=false;

    // load full results
    fetch('/api/dossier').then(r=>r.json()).then(data=>{
      const nAdv=data.summary?.n_advance||0;
      const nRej=data.summary?.n_reject||0;
      buildFunnel([
        {label:'측정 유전자',n:19944},
        {label:'차등발현 (FDR<0.01, |log2FC|>1)',n:3525},
        {label:'Cell-of-origin 종양세포 유래',n:100,cut:true},
        {label:'정상조직 안전성 통과',n:100-nRej,cut:true},
        {label:'ADVANCE',n:nAdv,final:true}
      ]);
      const allTargets=(data.advance_targets||[]).map(t=>({...t,judgment:'ADVANCE'}));
      buildCandidates(allTargets);
      buildRouting(allTargets);
      buildDossier(allTargets);
      if(data.contrasts)$('planSummary').innerHTML='<b>질환:</b> '+(data.disease||'TNBC')+' · <b>Contrasts:</b> '+data.contrasts.join(', ');
    }).catch(()=>{});

    // load omics charts
    fetch('/api/omics_charts').then(r=>r.json()).then(cd=>{
      if(cd.distance) drawDistance(cd.distance);
      if(cd.pca) drawPCA(cd.pca, cd.pca_var);
      if(cd.volcano) drawVolcano(cd.volcano);
      if(cd.gsea) drawGSEA(cd.gsea);
    }).catch(()=>{});
  }
}

async function startRun(){
  $('btnRun').disabled=true;
  $('statusChip').className='chip run';$('statusText').textContent='RUNNING';
  logN=0;evN=0;
  $('logBox').innerHTML='';$('funnelBox').innerHTML='<div class="locked">실행 중...</div>';
  $('candBox').innerHTML='<div class="locked">실행 중...</div>';
  $('routeBox').innerHTML='<div class="locked">실행 중...</div>';
  $('dosBox').innerHTML='<div class="locked">실행 중...</div>';
  ['tvDEG','tvCO','tvADV','tvREJ','tvEV','tvCR'].forEach(id=>$(id).textContent='—');
  AGENTS.forEach(a=>{$('nav-'+a.id).className='nv';$('st-'+a.id).textContent='IDLE';});
  STEPS.forEach((_,i)=>{const el=$('sp'+i);if(el)el.className='sp';});
  $('rvL').textContent='0';$('rvE').textContent='0';$('rvT').textContent='0s';
  ['riL','riE','riT'].forEach(id=>$(id).style.width='0');

  const st=Date.now();
  ti=setInterval(()=>{$('rvT').textContent=((Date.now()-st)/1000).toFixed(1)+'s';$('riT').style.width=Math.min(100,(Date.now()-st)/120000*100)+'%';},200);

  await fetch('/api/run',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({question:$('questionText').textContent})});
}

// ═══ 오믹스 차트 렌더링 ═══

function drawVolcano(data){
  $('rVol').textContent='완료';
  const cv=$('cvVol'),ctx=cv.getContext('2d');
  const W=cv.width,H=cv.height,pad={l:50,r:20,t:20,b:40};
  const pw=W-pad.l-pad.r,ph=H-pad.t-pad.b;
  ctx.clearRect(0,0,W,H);

  const xs=data.map(d=>d.x),ys=data.map(d=>d.y);
  const xMin=Math.min(...xs,-5),xMax=Math.max(...xs,5);
  const yMax=Math.min(Math.max(...ys),300);
  const sx=x=>(x-xMin)/(xMax-xMin)*pw+pad.l;
  const sy=y=>pad.t+ph-(y/yMax)*ph;

  // axes
  ctx.strokeStyle='#C4C6D0';ctx.lineWidth=1;
  ctx.beginPath();ctx.moveTo(pad.l,pad.t);ctx.lineTo(pad.l,H-pad.b);ctx.lineTo(W-pad.r,H-pad.b);ctx.stroke();
  // threshold lines
  ctx.strokeStyle='#E0E3E5';ctx.setLineDash([4,4]);
  ctx.beginPath();ctx.moveTo(sx(-1),pad.t);ctx.lineTo(sx(-1),H-pad.b);ctx.moveTo(sx(1),pad.t);ctx.lineTo(sx(1),H-pad.b);ctx.stroke();
  const y05=sy(-Math.log10(0.01));
  ctx.beginPath();ctx.moveTo(pad.l,y05);ctx.lineTo(W-pad.r,y05);ctx.stroke();
  ctx.setLineDash([]);

  // labels
  ctx.fillStyle='#757780';ctx.font='11px JetBrains Mono,monospace';ctx.textAlign='center';
  ctx.fillText('log2FC',W/2,H-5);
  ctx.save();ctx.translate(12,H/2);ctx.rotate(-Math.PI/2);ctx.fillText('-log10(FDR)',0,0);ctx.restore();

  // points
  data.forEach(d=>{
    const px=sx(d.x),py=sy(Math.min(d.y,yMax));
    ctx.beginPath();ctx.arc(px,py,2.5,0,Math.PI*2);
    ctx.fillStyle=d.s==='up'?'rgba(186,26,26,0.6)':d.s==='down'?'rgba(0,106,102,0.6)':'rgba(200,200,200,0.3)';
    ctx.fill();
  });

  // top gene labels
  const top=data.filter(d=>d.s!=='ns').sort((a,b)=>b.y-a.y).slice(0,6);
  ctx.font='9px Inter,sans-serif';ctx.fillStyle='#191C1E';ctx.textAlign='left';
  top.forEach(d=>{ctx.fillText(d.g,sx(d.x)+4,sy(Math.min(d.y,yMax))-4);});
}

function drawGSEA(data){
  $('rGsea').textContent='완료';
  const cv=$('cvGsea'),ctx=cv.getContext('2d');
  const W=cv.width,H=cv.height,pad={l:180,r:30,t:20,b:20};
  ctx.clearRect(0,0,W,H);

  const n=data.length;if(!n)return;
  const barH=Math.min(18,(H-pad.t-pad.b)/n-2);
  const maxNes=Math.max(...data.map(d=>Math.abs(d.nes)),2);

  data.forEach((d,i)=>{
    const y=pad.t+i*(barH+3);
    const bw=Math.abs(d.nes)/maxNes*(W-pad.l-pad.r);
    // bar
    ctx.fillStyle=d.nes>0?'rgba(0,106,102,0.7)':'rgba(186,26,26,0.5)';
    ctx.fillRect(pad.l,y,bw,barH);
    // NES label
    ctx.fillStyle='#191C1E';ctx.font='10px JetBrains Mono,monospace';ctx.textAlign='left';
    ctx.fillText(d.nes.toFixed(2),pad.l+bw+4,y+barH-3);
    // term label
    ctx.fillStyle='#44474F';ctx.font='10px Inter,sans-serif';ctx.textAlign='right';
    ctx.fillText(d.term,pad.l-6,y+barH-3);
  });
}

function drawDistance(data){
  $('rDist').textContent='완료';
  const cv=$('cvDist'),ctx=cv.getContext('2d');
  const W=cv.width,H=cv.height;
  ctx.clearRect(0,0,W,H);
  const n=data.n, mat=data.matrix, labels=data.labels;
  const pad={l:10,t:10,r:10,b:10};
  const cellW=(W-pad.l-pad.r)/n, cellH=(H-pad.t-pad.b)/n;

  // find max for color scale
  let mx=0;
  mat.forEach(row=>row.forEach(v=>{if(v>mx)mx=v;}));
  if(mx===0)mx=1;

  mat.forEach((row,i)=>{
    row.forEach((v,j)=>{
      const ratio=v/mx;
      // TNBC=warm, nonTNBC=cool
      const li=labels[i],lj=labels[j];
      const same=(li===lj);
      if(same){
        // low distance = dark teal
        const r2=Math.max(0,1-ratio);
        ctx.fillStyle=li==='TNBC'?`rgba(186,26,26,${0.1+r2*0.8})`:`rgba(0,106,102,${0.1+r2*0.8})`;
      } else {
        // cross-group = gray scale
        ctx.fillStyle=`rgba(117,119,128,${0.05+ratio*0.6})`;
      }
      ctx.fillRect(pad.l+j*cellW, pad.t+i*cellH, cellW, cellH);
    });
  });

  // border lines between TNBC / nonTNBC blocks
  const nTnbc=labels.filter(l=>l==='TNBC').length;
  ctx.strokeStyle='#191C1E';ctx.lineWidth=1.5;
  ctx.strokeRect(pad.l, pad.t, nTnbc*cellW, nTnbc*cellH);
  ctx.strokeRect(pad.l+nTnbc*cellW, pad.t+nTnbc*cellH, (n-nTnbc)*cellW, (n-nTnbc)*cellH);

  // labels
  ctx.fillStyle='#BA1A1A';ctx.font='bold 10px Inter';ctx.textAlign='left';
  ctx.fillText('TNBC',pad.l+2,pad.t+nTnbc*cellH+14);
  ctx.fillStyle='#006A66';
  ctx.fillText('non-TNBC',pad.l+nTnbc*cellW+2,H-pad.b-4);
}

function drawPCA(data, varExpl){
  $('rPca').textContent='완료';
  const cv=$('cvPca'),ctx=cv.getContext('2d');
  const W=cv.width,H=cv.height,pad={l:50,r:20,t:20,b:40};
  const pw=W-pad.l-pad.r,ph=H-pad.t-pad.b;
  ctx.clearRect(0,0,W,H);

  const xs=data.map(d=>d.x),ys=data.map(d=>d.y);
  const xMin=Math.min(...xs),xMax=Math.max(...xs);
  const yMin=Math.min(...ys),yMax=Math.max(...ys);
  const xRange=xMax-xMin||1, yRange=yMax-yMin||1;
  const sx=v=>pad.l+(v-xMin)/xRange*pw;
  const sy=v=>pad.t+ph-(v-yMin)/yRange*ph;

  // grid
  ctx.strokeStyle='#E6E8EA';ctx.lineWidth=0.5;
  for(let i=0;i<=4;i++){
    const gx=pad.l+pw*i/4;ctx.beginPath();ctx.moveTo(gx,pad.t);ctx.lineTo(gx,H-pad.b);ctx.stroke();
    const gy=pad.t+ph*i/4;ctx.beginPath();ctx.moveTo(pad.l,gy);ctx.lineTo(W-pad.r,gy);ctx.stroke();
  }

  // axes
  ctx.strokeStyle='#C4C6D0';ctx.lineWidth=1;
  ctx.beginPath();ctx.moveTo(pad.l,pad.t);ctx.lineTo(pad.l,H-pad.b);ctx.lineTo(W-pad.r,H-pad.b);ctx.stroke();
  const v1=varExpl?varExpl[0]:'?', v2=varExpl?varExpl[1]:'?';
  ctx.fillStyle='#757780';ctx.font='11px JetBrains Mono';ctx.textAlign='center';
  ctx.fillText(`PC1 (${v1}%)`,W/2,H-5);
  ctx.save();ctx.translate(14,H/2);ctx.rotate(-Math.PI/2);ctx.fillText(`PC2 (${v2}%)`,0,0);ctx.restore();

  // points
  data.forEach(d=>{
    const px=sx(d.x),py=sy(d.y);
    ctx.beginPath();ctx.arc(px,py,4,0,Math.PI*2);
    ctx.fillStyle=d.g==='TNBC'?'rgba(186,26,26,0.65)':'rgba(0,106,102,0.5)';
    ctx.fill();
    ctx.strokeStyle='#fff';ctx.lineWidth=0.5;ctx.stroke();
  });

  // legend
  [{c:'#BA1A1A',l:'TNBC'},{c:'#006A66',l:'non-TNBC'}].forEach((lg,i)=>{
    ctx.fillStyle=lg.c;ctx.beginPath();ctx.arc(W-90,pad.t+i*20+10,5,0,Math.PI*2);ctx.fill();
    ctx.fillStyle='#44474F';ctx.font='11px Inter';ctx.textAlign='left';ctx.fillText(lg.l,W-80,pad.t+i*20+14);
  });
}

conn();
</script>
</body>
</html>"""
