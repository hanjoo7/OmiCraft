/* Recorded pipeline report. No synthetic measurements or inference requests. */
function renderPipelineScenario(d) {
 'use strict';
 const $=id=>document.getElementById(id), el=(p,t,v,c)=>{const n=document.createElement(t);if(v!=null)n.textContent=String(v);if(c)n.className=c;p.append(n);return n;};
 const val=v=>v==null?'미기록':typeof v==='object'?JSON.stringify(v):String(v);
 const count=v=>v==null?'—':Number(v).toLocaleString('ko-KR');
 const s=d.summary||{}, details=d.pipeline_details||{}, plan=details.plan||{}, qual=details.qualification||{};
 const hiddenGenes=new Set(d.candidate_display?.hidden_genes||[]),visible=c=>!hiddenGenes.has(String(c.gene||'').trim().toUpperCase());
 const candidates=(s.candidates||[]).filter(visible), evidence=details.evidence||[], routes=(s.modalities||[]).filter(visible), design=d.design_report||{};
 const stages=(d.stages||[]).filter(r=>r.stage!=='readiness'), budget=plan.budget||{};
 const cls=v=>['ADVANCE','PASS','COMPLETED','APPROVE'].includes(v)?'adv':['REJECT','FAIL','FAILED'].includes(v)?'rej':'hold';
 const box=(p,title)=>{const b=el(p,'div',null,'card');if(title)el(b,'div',title,'ch');return b;};
 const text=(p,t,c)=>el(p,'p',t,c||'pr-note');
 const link=(p,label,url)=>{const a=el(p,'a',label);if(url&&(/^(?:\/[^/]|\.\.?\/|https?:\/\/)/i.test(url)))a.href=url;return a;};
 const raw=(p,label,data)=>{const x=el(p,'details',null,'pr-raw');el(x,'summary',label);el(x,'pre',JSON.stringify(data,null,2));return x;};
 const table=(p,rows,fields)=>{const wrap=el(p,'div',null,'pr-table');const t=el(wrap,'table'),h=el(el(t,'thead'),'tr'),b=el(t,'tbody');fields.forEach(f=>el(h,'th',f[0]));rows.forEach(r=>{const tr=el(b,'tr');fields.forEach(f=>el(tr,'td',val(typeof f[1]==='function'?f[1](r):r[f[1]])));});return wrap;};
 document.title='OmiCraft — Pipeline Report · '+d.run_id;
 renderReportNavigation(design.reference,'pipeline');
 const successfulGenes=new Set([...(design.panels||[]),...(design.followups||[])].filter(designRoutePassed).map(p=>p.gene));
 // 01 — input and constraints
 const q=el($('pr-question'),'div',null,'qbox'), qb=box(q,'RESEARCH QUESTION'), body=el(qb,'div',null,'qbody');
 el(body,'div',s.question||d.question||'연구질문 미기록','qq');
 text(body,(plan.cohort_description||'코호트 설명 미기록')+' · '+(s.de?.contrast||'분석 contrast 미기록'));
 const manifest=s.input_manifest||{}, constraints=box(q,'INPUT / CONSTRAINT'), dl=el(constraints,'dl',null,'qdl');
 for(const [k,v] of [['분석 샘플',count(manifest.samples??s.survival?.samples)],['분석군',Object.entries(manifest.group_counts||{}).map(([k,v])=>k+' '+count(v)).join(' · ')||'미기록'],['후보 예산',budget.max_candidates],['GPU 예산',budget.gpu_hours==null?null:budget.gpu_hours+'h'],['설계 선택',s.selection?.gate_status]]){el(dl,'dt',k);el(dl,'dd',val(v));}
 raw($('pr-question'),'실행 계획 · 비교군 · 중단 조건',{plan,recorded_contrasts:s.contrasts,input_manifest:manifest});
 // Recorded resource accounting in the sidebar
 function resource(label,value,max){const row=el($('pr-resources'),'div',null,'srb'),l=el(row,'div',null,'l');el(l,'span',label);el(l,'b',value==null?'미기록':val(value)+(max?' / '+max:''));const bar=el(el(row,'div',null,'t'),'i');bar.style.width=value!=null&&max?Math.min(100,value/max*100)+'%':'0%';}
 const usage=s.resource_usage||d.resource_usage||{};
 resource('LLM Tokens',usage.total_tokens,budget.llm_tokens??budget.llm_token_budget);resource('GPU 할당 시간 (h)',usage.gpu_hours==null?null:Number(usage.gpu_hours.toFixed(4)),budget.gpu_hours);resource('계측 도구 호출',usage.tool_calls,budget.tool_calls);
 if(usage.schema_version){resource('결과 재사용',usage.cache_hits);text($('pr-resources'),'이번 실행의 API 토큰·분석/설계 도구 호출·GPU 할당 시간을 기록합니다. GPU 대기 시간과 이전 실행의 계산량은 제외합니다.','pr-muted');if(usage.llm_usage_missing)text($('pr-resources'),'확인된 토큰 '+count(usage.recorded_tokens)+' · 사용량 응답 미확인 '+usage.llm_usage_missing+'건','pr-muted');if(usage.unmetered_children?.length)text($('pr-resources'),'자원 기록이 없는 설계 실행 '+usage.unmetered_children.length+'건 · 전체 합계 미확정','pr-muted');if(usage.gpu_unfinished_leases)text($('pr-resources'),'종료되지 않은 GPU 할당 '+usage.gpu_unfinished_leases+'건 · 종료 후 합계 확정','pr-muted');}else text($('pr-resources'),'자원 계측 도입 전 실행입니다. 기록이 없는 과거 사용량은 복원할 수 없습니다.','pr-muted');
 if(s.r_execution_mode)text($('pr-resources'),'R 분석: '+(s.r_execution_mode==='VERIFIED_CACHE_REUSE'?'검증된 완료 결과 재사용':s.r_execution_mode==='FRESH_ANALYSIS'?'새로 계산':s.r_execution_mode),'pr-muted');
 const elapsed=stages.length>1?(Date.parse(stages.at(-1).timestamp)-Date.parse(stages[0].timestamp))/1000:null;
 if(Number.isFinite(elapsed)&&elapsed>=0)text($('pr-resources'),'기록 시간 간격 '+Math.round(elapsed)+'초','pr-muted');
 $('pr-status').className='chip '+(d.running?'run':d.execution_status==='COMPLETED'?'ok':['FAILED','BLOCKED'].includes(d.execution_status)?'failed':'warn');
 $('pr-status').textContent=d.running?'실행 중 · 기록 갱신':'실행 '+d.execution_status;
 $('pr-print').onclick=()=>window.print();
 // 02 — recorded figures after batch correction
 const plots=(d.artifacts||[]).filter(a=>a.artifact_type==='plot'&&!/(^|\/)before_batch_/i.test(a.path)), used=new Set();
 function figure(parent,a,label){const b=box(parent,label||a?.title||'산출물 미기록');if(!a){text(b,'시각화 불가: 이 실행에 해당 그림 파일이 없습니다.','pr-empty');return;}used.add(a.path);if(!a.available||!a.href){text(b,'시각화 불가: 저장된 파일을 찾을 수 없거나 접근할 수 없습니다.','pr-empty');return;}const img=el(b,'img');img.src=a.href;img.alt=a.title||a.label||label;img.loading='lazy';img.className='pr-plot';img.onerror=()=>{img.hidden=true;text(b,'시각화 불가: 이미지를 불러오지 못했습니다.','pr-empty');};text(b,a.caption||a.title||'', 'cnote');}
 const op=$('pr-omics');
 const corrected=el(op,'div',null,'charts');for(const [suffix,title] of [['sample_distance.png','Sample distance matrix'],['pca_group.png','PCA · Tumor group'],['pca_tss.png','PCA · Tissue source site']])figure(corrected,plots.find(a=>a.path.endsWith('after_batch_'+suffix)),title+' · 배치 보정 후');
 const charts=el(op,'div',null,'charts');if(details.volcano?.points?.length){renderRecordedVolcano(box(charts,'VOLCANO — DIFFERENTIAL EXPRESSION'),details.volcano,s.de||{});plots.filter(a=>a.path.endsWith('/volcano.png')).forEach(a=>used.add(a.path));}else figure(charts,plots.find(a=>a.path.endsWith('/volcano.png')),'VOLCANO — DIFFERENTIAL EXPRESSION');figure(charts,plots.find(a=>a.path.endsWith('/Hallmark_NES_top.png')),'GSEA — HALLMARK');
 text(op,'GSEA는 저장된 NES 요약 그림을 표시합니다. Running enrichment curve가 저장되지 않은 실행에서는 곡선을 생성하지 않습니다.');
 const networks=plots.filter(a=>/GOBP_network[^/]*\.png$/i.test(a.path)||a.title==='GO Biological Process Network');
 const network=el(op,'div',null,'pr-network');for(const a of networks)figure(network,a,'GO Biological Process Network');
 const other=el(op,'details');el(other,'summary','추가 분석 그림 · 경로 농축 · 발현 · 세포 맥락 · dependency');const rest=el(other,'div',null,'charts');plots.filter(a=>!used.has(a.path)).forEach(a=>figure(rest,a));
 raw(op,'DE · GSEA · Survival 통계',{de:s.de,gsea:s.gsea,survival:s.survival});
 // 03 — comparable observed counts, no fabricated external-validation funnel
 const directions=s.de?.direction_counts||{},deg=Object.keys(directions).length?(directions.Up_in_TNBC||0)+(directions.Down_in_TNBC||0):null;
 const funnel=[['측정 유전자',manifest.genes_raw,'입력 manifest'],['발현 필터 통과',manifest.genes_de_and_survival,'DE / 생존 분석 유전자'],['차등발현 유전자',deg,'DESeq2 · FDR '+val(s.de?.alpha)+' · |log2FC| '+val(s.de?.absolute_log2FC)],['적격성 평가 대상',candidates.length||null,'설정된 후보 예산 내 평가'],['ADVANCE 표적',candidates.length?candidates.filter(c=>c.decision==='ADVANCE').length:null,'표적 수준의 판정'],['계산 완료 표적',routes.length?new Set(routes.filter(r=>r.execution_status==='COMPLETED').map(r=>r.gene)).size:null,'원래 실행의 고비용 설계'],['분자 검증 통과 경로',design.passing_routes??null,'원래 실행 · 후속 결과는 표적별 설계 결과에서 확인']];
 const fb=el(box($('pr-funnel'),'CANDIDATE FUNNEL · RECORDED COUNTS'),'div',null,'funnel'),max=Math.max(1,...funnel.map(r=>r[1]||0));
 for(const [label,n,note] of funnel){const row=el(fb,'div',null,'fr on');const l=el(row,'div',label,'fl');el(l,'span',note,'drop');const bar=el(el(row,'div',null,'ft'),'i');bar.style.width=n==null?'0%':Math.log1p(n)/Math.log1p(max)*100+'%';el(row,'div',count(n),'fn');}
 text($('pr-funnel'),'막대 길이는 log(1+n) 스케일입니다. 각 단계의 단위와 선정 기준이 다르며, 후속 검증은 원래 실행의 통과 수에 합산하지 않습니다.');
 // 04 — actionable decisions; HOLD records remain in target details and raw data
 const qp=$('pr-qualification'),tools=el(qp,'div',null,'pr-filters'),search=el(tools,'input');search.type='search';search.placeholder='표적명 검색';search.setAttribute('aria-label','표적명 검색');search.id='pr-search';
 const filter=el(tools,'select');filter.setAttribute('aria-label','판정 필터');for(const v of ['ALL','ADVANCE','REJECT']){const o=el(filter,'option',v);o.value=v;}
 const totals=el(tools,'span',null,'pr-muted'),cards=el(qp,'div',null,'cands');
 const qualificationCandidates=candidates.filter(c=>c.decision!=='HOLD'),qualificationGenes=new Set(qualificationCandidates.map(c=>c.gene));
 const qualificationEvidence=evidence.filter(e=>qualificationGenes.has(e.gene));
 const evDetail=el(qp,'details');el(evDetail,'summary','Evidence Store · '+qualificationEvidence.length+'건');const evGrid=el(evDetail,'div',null,'evgrid'),more=el(evDetail,'button','근거 더 보기','btn ghost');let evLimit=24;
 let selectedGene=null;
 const sorted=[...candidates].sort((a,b)=>Number(successfulGenes.has(b.gene))-Number(successfulGenes.has(a.gene))||(a.decision==='ADVANCE'?0:1)-(b.decision==='ADVANCE'?0:1));
 function matches(c){return (!search.value||String(c.gene).toLowerCase().includes(search.value.toLowerCase()))&&(filter.value==='ALL'||c.decision===filter.value);}
 function paintCandidates(){
  cards.replaceChildren();const list=sorted.filter(c=>c.decision!=='HOLD').filter(matches);totals.textContent=list.length+' / '+qualificationCandidates.length+'개 표적';
  for(const c of list){const b=el(cards,'button',null,'cd live '+cls(c.decision));b.type='button';b.dataset.gene=c.gene;b.setAttribute('aria-pressed',String(c.gene===selectedGene));b.onclick=()=>{selectTarget(c.gene);$('s-dossiers').scrollIntoView({block:'start'});};
   const head=el(b,'div',null,'cht');el(head,'span',c.gene,'cg');el(head,'span',c.decision,'cv');if(successfulGenes.has(c.gene)){b.classList.add('pr-success-target');el(b,'span','설계 계산 검증 통과','ir-result-badge');}el(b,'div','Biological '+val(c.B)+' · Tier '+val(c.tier),'cs');const grade=el(b,'div',null,'grade');grade.setAttribute('aria-label','Biological confidence '+val(c.B));const level={B1:3,B2:2,B3:1}[c.B]||0;for(let i=1;i<=3;i++)el(grade,'i',null,i<=level?'on':'');el(b,'div',c.reason,'cw');el(b,'div',(c.routes||[]).map(r=>r.modality+' '+r.Tier).join(' · '),'cm');}
  if(!list.length)text(cards,'조건에 맞는 표적 기록이 없습니다.','pr-empty');
  evGrid.replaceChildren();const genes=new Set(list.map(c=>c.gene)),ev=qualificationEvidence.filter(e=>genes.has(e.gene));
  for(const e of ev.slice(0,evLimit)){const b=el(evGrid,'div',null,'ev '+({SUPPORTIVE:'sup',CONTRADICTORY:'con',NOT_APPLICABLE:'na'}[e.status]||'neu'));const h=el(b,'div',null,'eh');el(h,'span',e.gene,'eid');el(h,'span',e.status,'ety');el(b,'div',e.claim,'ec');el(b,'div',e.source+' · '+e.id,'es');raw(b,'근거 원문 · 출처',e);}
  more.hidden=ev.length<=evLimit;more.textContent='근거 더 보기 ('+Math.min(evLimit,ev.length)+' / '+ev.length+')';
 }
 search.oninput=()=>{evLimit=24;paintCandidates();};filter.onchange=()=>{evLimit=24;paintCandidates();};more.onclick=()=>{evLimit+=48;paintCandidates();};paintCandidates();
 // 05 — actual design routes and recorded reasons
 const router=box($('pr-routing'),'MODALITY ROUTER · 4 BRANCHES');
 function routeReason(r){const b=r.blocker||'';if(b.includes('target_safety_veto'))return '안전성 기준 미충족';if(b.includes('BLOCKED_TIER3'))return 'Tier 3 · 자동 설계 진입 제한';if(b.includes('BLOCKED_HOLD'))return '필수 근거 부족 · HOLD';if(b.includes('antibody_antigen_complex'))return '항체–항원 복합체 및 항원 chain 입력 필요';if(b.includes('warhead_e3_linker'))return 'Warhead·E3·linker 및 연결 위치 입력 필요';return b||(r.execution_status==='COMPLETED'?'계산 완료':'사유 미기록');}
 const branches=[['Antibody / ADC',['ADC'],'세포 표면 접근 · 항체/항원 복합체'],['Small molecule',['SMALL_MOLECULE'],'세포 내 효소 · pocket · ligand'],['Degrader',['DEGRADER'],'target binder · E3 · linker'],['Binder / ligand trap',['DE_NOVO_BINDER'],'설계 가능한 표적 구조 · 결합 위치']];
 for(const [name,mods,why] of branches){const rr=routes.filter(r=>mods.includes(r.modality)),done=rr.filter(r=>r.execution_status==='COMPLETED').length;const row=el(router,'details',null,'pr-route');const h=el(row,'summary');el(h,'span',name);el(h,'span',!mods.length?'계산 경로 미구현':!rr.length?'실행 기록 없음':done+' 완료 / '+rr.length+' 경로','pr-muted');text(row,why);const reasons={};rr.filter(r=>r.execution_status!=='COMPLETED').forEach(r=>{const reason=routeReason(r);reasons[reason]=(reasons[reason]||0)+1;});if(Object.keys(reasons).length)text(row,Object.entries(reasons).map(([reason,n])=>reason+' '+n+'건').join(' · '),'pr-route-reasons');if(rr.length)table(row,rr,[['표적','gene'],['실행','execution_status'],['검증','validation_status'],['Critic','critic_decision'],['실행 방식','execution_mode'],['미실행 / 차단 사유',routeReason]]);else text(row,!mods.length?'분기는 표시되지만 이 버전에는 해당 계산 workflow가 구현되어 있지 않습니다.':'이 실행에서 선택·실행된 경로가 없습니다.');}
 text($('pr-routing'),'구현된 4개 모달리티의 실제 실행·차단 사유를 각 분기에서 확인할 수 있습니다. De novo binder는 중화 기능이 확인되었다는 의미가 아닙니다.');
 // 06 — selected target dossier and suggested, explicitly unexecuted wet-lab plans
 const wet={ADC:['항원 결합 및 내재화 평가','정상조직 counter-screen','접합체 안정성·DAR·세포독성 평가'],DE_NOVO_BINDER:['발현·정제 및 응집 평가','SPR/BLI 결합과 비표적 결합 평가','세포 기반 기능 및 선택성 검증'],SMALL_MOLECULE:['생화학적 활성·결합 평가','세포 내 target engagement 및 선택성 평가','노출·독성 및 ADMET 검증'],DEGRADER:['표적 분해 DC50/Dmax 및 기전 확인','E3·proteasome 의존성 확인','Proteomics 선택성 평가']};
 const dos=$('pr-dossiers'),pickerLabel=el(dos,'label','표적 선택 ','pr-target-picker'),picker=el(pickerLabel,'select');picker.id='pr-target';picker.setAttribute('aria-label','결과를 볼 표적');for(const c of sorted){const o=el(picker,'option',c.gene+' · '+c.decision+(successfulGenes.has(c.gene)?' · ★ 설계 검증 통과':''));o.value=c.gene;}if(design.reference){const o=el(picker,'option','ERBB2 · 이전 실행 참고');o.value='__reference__';}const success=el(dos,'div',null,'pr-success-targets');
 for(const c of sorted.filter(c=>successfulGenes.has(c.gene))){const button=el(success,'button','★ '+c.gene+' · 설계 검증 통과','btn ghost');button.type='button';button.onclick=()=>selectTarget(c.gene);}
 if(success.children.length)text(success,'Validation·Critic이 모두 ADVANCE인 완료 결과입니다. 후속 결과를 포함하며 실험 검증과는 별도입니다.','pr-muted');
 const targetBody=el(dos,'div');text(dos,'실험 계획은 모달리티별 후속 검증 제안이며, 실행된 실험이나 검증 결과가 아닙니다.');
 function selectTarget(gene){
  selectedGene=gene;picker.value=gene;targetBody.replaceChildren();$('pr-design').replaceChildren();document.querySelectorAll('#pr-qualification .cd').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.gene===gene)));
  const c=sorted.find(c=>c.gene===gene);if(c){const info=qual[c.gene]||{},b=el(targetBody,'div',null,'card dossier');b.id='dossier-'+encodeURIComponent(c.gene);const h=el(b,'div',null,'dh');el(h,'b',c.gene);el(h,'span',c.decision+' · '+val(c.tier),'dtag');const rows=el(b,'div',null,'drows');
  for(const [k,v] of [['BIOLOGICAL CONFIDENCE',c.B],['INTERVENTION',info.intervention_role],['LOCATION',info.target_info?.location],['MISSING EVIDENCE',info.missing_evidence],['REASON',c.reason]]){const r=el(rows,'div',null,'drow');el(r,'div',k,'dk');el(r,'div',val(v),'dv');}
  table(b,c.routes||[],[['모달리티','modality'],['B','B'],['M','M'],['D','D'],['Tier','Tier'],['적합성 근거','modality_reason'],['개발 근거','readiness_reason'],['Safety veto',r=>r.safety_veto?.has_veto?'REJECT · '+r.safety_veto.veto_reason:(r.safety_veto?.hold_reasons||[]).join(' · ')||'기록된 veto 없음']]);
  const ev=evidence.filter(e=>e.gene===c.gene);raw(b,'지지·반대 근거 '+ev.length+'건',ev);
  const w=el(b,'div',null,'dwet');el(w,'div','WET-LAB · 후속 검증 제안 / NOT EVALUATED','dk');const mods=[...new Set(routes.filter(r=>r.gene===c.gene&&r.execution_status==='COMPLETED').map(r=>r.modality))];if(!mods.length){
   const recorded=routes.filter(r=>r.gene===c.gene),active=(design.panels||[]).some(p=>p.gene===c.gene&&(p.running||p.execution_status==='RUNNING'))||recorded.some(r=>r.execution_status==='RUNNING');
   if(active)text(w,c.gene+' 설계가 진행 중입니다. 완료 후 해당 경로의 후속 실험 계획을 표시합니다.');
   else if(recorded.length){
    text(w,c.gene+'에서 계산이 완료된 설계 경로가 없어 후속 실험 계획을 표시하지 않습니다.');
    for(const r of recorded)text(w,r.modality+' · '+(r.execution_status==='BLOCKED'?'설계 진입 차단':r.execution_status==='FAILED'?'실행 실패':val(r.execution_status))+' · '+routeReason(r));
   }else text(w,c.gene+'의 설계 실행 기록이 아직 없습니다. 표적 판정: '+val(c.decision)+'.');
  }for(const m of mods){el(w,'b',m);const ol=el(w,'ol');(wet[m]||[]).forEach(t=>el(ol,'li',t));}

 }
  if(gene==='__reference__')text(targetBody,'ERBB2 · 이전 실행의 결과를 아래에서 확인할 수 있습니다.');
  const panels=(design.panels||[]).filter(p=>p.gene===gene),followups=(design.followups||[]).filter(p=>p.gene===gene);
  const completed=panels.filter(p=>!p.running&&p.execution_status==='COMPLETED');
  if(d.design_report)renderPipelineDesign($('pr-design'),{...design,panels,followups,candidate_search:{},evaluated_targets:c?1:0,candidate_count:panels.reduce((n,p)=>n+(p.candidate_count||0),0),completed_routes:completed.length,passing_routes:completed.filter(p=>p.validation==='ADVANCE'&&p.critic==='ADVANCE'&&!p.is_mock).length,reference:gene==='__reference__'?design.reference:null},{downloads:false,gene:c?.gene,referenceOnly:gene==='__reference__'});
  else text($('pr-design'),'저장된 설계 결과가 없습니다.');
 }
 picker.onchange=()=>selectTarget(picker.value);
 if(sorted.length)selectTarget(sorted[0].gene);else if(design.reference)selectTarget('__reference__');else text(targetBody,'표적 결과를 생성할 적격성 기록이 없습니다.','pr-empty');
 // 08 — provenance and recorded review
 const files=$('pr-files');text(files,'실행 ID · '+d.run_id);
 for(const msg of [...(s.limitations||[]),...(d.report_warnings||[])])text(files,msg,'pr-warning');
 raw(files,'Agent interpretations · 과학적 검토 의견', {notes:s.agent_notes,issues:s.review_issues,blockers:d.blockers,errors:d.errors});
 raw(files,'전체 실행 기록 · 통계 · 판정 원본',{...d,artifacts:d.artifacts});
 const reused=el(files,'div');for(const row of s.reused_stages||[]){text(reused,'분석 단계 재사용 · '+row.stage);link(reused,row.source_run_id,'/report?run_id='+encodeURIComponent(row.source_run_id));}
 if(d.running&&!window.__COMPETITION_RUN__){window.clearTimeout(window.omicraftReportTimer);window.omicraftReportTimer=setTimeout(()=>location.reload(),10000);}
}

/* Full recorded observations, with a selectable viewport rather than outlier-driven axes. */
function renderRecordedVolcano(parent,data,de){
 const el=(p,t,v)=>{const n=document.createElement(t);if(v!=null)n.textContent=v;p.append(n);return n;};
 const controls=el(parent,'div');controls.className='pr-filters';
 function selector(label,id,values,initial){const l=el(controls,'label',label+' '),s=el(l,'select');s.id=id;for(const [v,t] of values){const o=el(s,'option',t);o.value=v;}s.value=initial;return s;}
 const xs=selector('log2FC 범위','pr-volcano-x',[['2','−2 ~ 2'],['4','−4 ~ 4'],['8','−8 ~ 8'],['all','전체']],'4');
 const ys=selector('−log10(padj)','pr-volcano-y',[['50','0 ~ 50'],['100','0 ~ 100'],['all','전체']],'100');
 const all=el(controls,'button','전체 범위');all.type='button';all.className='btn ghost';all.id='pr-volcano-all';
 const canvas=el(parent,'canvas');canvas.id='pr-volcano';canvas.width=1440;canvas.height=960;canvas.style.width='100%';canvas.style.height='auto';canvas.setAttribute('role','img');canvas.setAttribute('aria-label','실제 차등발현 유전자 volcano plot. 축 범위는 위 선택기로 조정할 수 있습니다.');
 const note=el(parent,'p');note.className='cnote';note.id='pr-volcano-count';note.setAttribute('aria-live','polite');
 const hover=el(parent,'p','점에 마우스를 올리면 유전자와 값을 확인할 수 있습니다.');hover.className='cnote';
 const legend=el(parent,'p','● 상향 (teal) · ● 하향 (navy) · ● 비유의 (회색)');legend.className='cnote';
 const points=data.points,ctx=canvas.getContext('2d'),W=720,H=480,L=76,R=694,T=24,B=410;
 const full=points.reduce((a,p)=>[Math.min(a[0],p[1]),Math.max(a[1],p[1]),Math.max(a[2],p[2])],[0,0,0]);
 const colors={Up_in_TNBC:'#006A66',Down_in_TNBC:'#3E5688',Not_DEG:'#B4B9BE'};let shown=[];
 function draw(){
  const xmin=xs.value==='all'?Math.floor(full[0])-1:-Number(xs.value),xmax=xs.value==='all'?Math.ceil(full[1])+1:Number(xs.value),ymax=ys.value==='all'?Math.max(10,Math.ceil(full[2]/10)*10):Number(ys.value);
  const px=x=>L+(x-xmin)/(xmax-xmin)*(R-L),py=y=>B-y/ymax*(B-T);
  ctx.setTransform(2,0,0,2,0,0);ctx.clearRect(0,0,W,H);ctx.fillStyle='#fff';ctx.fillRect(0,0,W,H);ctx.font='12px sans-serif';
  for(let i=0;i<=4;i++){const x=L+i*(R-L)/4,y=T+i*(B-T)/4;ctx.strokeStyle='#ECEEF0';ctx.beginPath();ctx.moveTo(x,T);ctx.lineTo(x,B);ctx.moveTo(L,y);ctx.lineTo(R,y);ctx.stroke();ctx.fillStyle='#757780';ctx.textAlign='center';ctx.fillText((xmin+i*(xmax-xmin)/4).toFixed(1),x,B+20);ctx.textAlign='right';ctx.fillText((ymax-i*ymax/4).toFixed(0),L-9,y+4);}
  ctx.save();ctx.beginPath();ctx.rect(L,T,R-L,B-T);ctx.clip();shown=[];
  const groups={Not_DEG:[],Down_in_TNBC:[],Up_in_TNBC:[]};let outside=0;
  for(const p of points){if(p[1]<xmin||p[1]>xmax||p[2]>ymax){outside++;continue;}const dot=[px(p[1]),py(p[2]),p];shown.push(dot);(groups[p[3]]||groups.Not_DEG).push(dot);}
  for(const [direction,dots] of Object.entries(groups)){ctx.fillStyle=colors[direction];ctx.globalAlpha=direction==='Not_DEG'?.18:.45;ctx.beginPath();for(const [x,y] of dots){ctx.moveTo(x+1.4,y);ctx.arc(x,y,1.4,0,2*Math.PI);}ctx.fill();}
  ctx.globalAlpha=1;ctx.strokeStyle='#9A5B06';ctx.setLineDash([5,4]);ctx.beginPath();const threshold=Number(de.absolute_log2FC),alpha=Number(de.alpha);if(Number.isFinite(threshold)&&threshold>0)for(const x of [-threshold,threshold]){ctx.moveTo(px(x),T);ctx.lineTo(px(x),B);}if(alpha>0&&alpha<1){const y=py(-Math.log10(alpha));ctx.moveTo(L,y);ctx.lineTo(R,y);}ctx.stroke();ctx.restore();
  ctx.strokeStyle='#C4C6D0';ctx.strokeRect(L,T,R-L,B-T);ctx.fillStyle='#44474F';ctx.textAlign='center';ctx.font='14px sans-serif';ctx.fillText('log2 fold change', (L+R)/2,H-20);ctx.save();ctx.translate(23,(T+B)/2);ctx.rotate(-Math.PI/2);ctx.fillText('−log10 adjusted p-value',0,0);ctx.restore();
  note.textContent=shown.length.toLocaleString()+' / '+points.length.toLocaleString()+'개 표시 · 범위 밖 '+outside.toLocaleString()+'개'+(data.skipped?' · 좌표 미기록 '+data.skipped+'개 제외':'')+'. 원본 값은 유지하며 전체 범위에서 모든 유효 점을 확인할 수 있습니다.';
  canvas.dataset.visible=shown.length;canvas.dataset.total=points.length;
 }
 xs.onchange=ys.onchange=draw;all.onclick=()=>{xs.value=ys.value='all';draw();};
 canvas.onmousemove=e=>{const r=canvas.getBoundingClientRect(),x=(e.clientX-r.left)*W/r.width,y=(e.clientY-r.top)*H/r.height;let nearest=null,distance=64;for(const dot of shown){const v=(dot[0]-x)**2+(dot[1]-y)**2;if(v<distance){nearest=dot[2];distance=v;}}hover.textContent=nearest?(nearest[0]||'이름 미기록')+' · log2FC '+nearest[1].toFixed(3)+' · −log10(padj) '+nearest[2].toFixed(3):'점에 마우스를 올리면 유전자와 값을 확인할 수 있습니다.';};
 draw();
}
