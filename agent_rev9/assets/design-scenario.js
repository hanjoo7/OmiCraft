/* Fixed-target design reports use the same shell as the discovery pipeline. */
function renderDesignScenario(d){
 'use strict';
 const $=id=>document.getElementById(id);
 const el=(p,t,v,c)=>{const n=document.createElement(t);if(v!=null)n.textContent=String(v);if(c)n.className=c;p.append(n);return n;};
 const value=v=>v==null?'미기록':typeof v==='object'?JSON.stringify(v):String(v);
 const note=(p,v)=>el(p,'p',v,'pr-note');
 const box=(p,title)=>{const b=el(p,'div',null,'card');el(b,'div',title,'ch');return b;};
 const raw=(p,title,data)=>{const b=el(p,'details',null,'pr-raw');el(b,'summary',title);el(b,'pre',JSON.stringify(data,null,2));};
 const report=d.integrated_report||{},panels=report.panels||[],gene=d.gene||'표적 미기록';
 const passed=designRoutePassed;
 renderReportNavigation({report_url:'/report?run_id='+encodeURIComponent(d.run_id)},'design');
 const tone=p=>passed(p)?'adv':['REJECT','FAIL','FAILED'].includes(p.validation)||p.execution_status==='FAILED'?'rej':'hold';
 const mode=p=>p.execution_mode==='REUSED_RESULT'?'기존 결과 재사용':p.execution_mode==='FRESH_ANALYSIS'?'새로 계산':value(p.execution_mode);
 document.title='OmiCraft — '+gene+' · Therapeutic Design Report';
 const status=$('pr-status');status.textContent=d.running?'실행 중 · 기록 갱신':'실행 '+value(d.execution_status);
 status.className='chip '+(d.running?'run':d.execution_status==='COMPLETED'?'ok':'warn');
 // 01: only recorded design inputs; no inferred discovery scores or cohorts.
 const q=el($('pr-question'),'div',null,'qbox'),body=el(box(q,'THERAPEUTIC DESIGN'),'div',null,'qbody');
 el(body,'h1',gene+' · 모달리티 설계 결과','qq');
 note(body,'Small molecule · De novo binder · ADC · Degrader');
 const dl=el(box(q,'INPUT / EXECUTION'),'dl',null,'qdl');
 for(const [k,v] of [['표적',gene],['실행 유형',d.execution_profile],['모달리티',d.modality],['실행 ID',d.run_id]]){el(dl,'dt',k);el(dl,'dd',value(v));}
 // Missing resource accounting stays missing rather than becoming zero.
 const usage=d.resource_usage||d.summary?.resource_usage||{};
 for(const [label,v] of [['LLM Tokens',usage.total_tokens],['GPU 할당 시간 (h)',usage.gpu_hours],['계측 도구 호출',usage.tool_calls]]){
  const row=el($('pr-resources'),'div',null,'srb'),line=el(row,'div',null,'l');el(line,'span',label);el(line,'b',value(v));
 }
 if(!usage.schema_version)note($('pr-resources'),'이 실행에는 자원 사용량 계측 기록이 없습니다.');
 // 02: execution and scientific decisions remain separate.
 const tiles=el($('pr-overview'),'div',null,'tiles');
 for(const [label,v] of [['설계 경로',panels.length],['계산 완료',panels.filter(p=>p.execution_status==='COMPLETED'&&!p.running).length],['계산 기준 통과',panels.filter(passed).length],['HOLD 판정',panels.filter(p=>p.validation==='HOLD').length]]){
  const t=el(tiles,'div',null,'tile');el(t,'div',label,'tl');el(t,'div',v,'tv');
 }
 note($('pr-overview'),'계산 기준 통과는 Validation과 Critic이 모두 ADVANCE인 실제 계산 경로입니다. 추가 품질 필터와 실험 검증은 상세 결과에서 별도로 확인합니다.');
 const cards=el($('pr-overview'),'div',null,'cands');
 const picker=el($('pr-design'),'nav',null,'pr-modality-picker');picker.id='pr-modality';picker.setAttribute('aria-label','결과를 볼 모달리티');
 const host=el($('pr-design'),'div',null,'integrated-report ir-compact-figures'),sections=[];
 function select(name){
  picker.dataset.selected=name;
  picker.querySelectorAll('button').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.modality===name)));
  for(const [key,section] of sections)section.hidden=key!==name;
  cards.querySelectorAll('button').forEach(b=>b.setAttribute('aria-pressed',String(b.dataset.modality===name)));
 }
 for(const p of panels){
  const card=el(cards,'button',null,'cd live '+tone(p));card.type='button';card.dataset.modality=p.modality;card.setAttribute('aria-pressed','false');
  const h=el(card,'div',null,'cht');el(h,'span',p.title,'cg');el(h,'span',p.validation||'미기록','cv');if(passed(p)){card.classList.add('pr-success-target');el(card,'span','★ 계산 검증 통과','ir-result-badge');}
  el(card,'div','실행 '+value(p.execution_status)+' · Critic '+value(p.critic),'cs');
  el(card,'div',p.scope,'cw');el(card,'div',mode(p),'cm');
  card.onclick=()=>{select(p.modality);$('s-design').scrollIntoView({block:'start'});};
  const button=el(picker,'button',p.title+(passed(p)?' · ★':''),'btn ghost');button.type='button';button.dataset.modality=p.modality;button.onclick=()=>select(p.modality);
  const section=el(host,'section',null,'ir-section');section.id='result-'+p.modality;sections.push([p.modality,section]);
  renderDesignPanel(section,p,{downloads:false});
 }
 if(!panels.length)note(host,'저장된 모달리티 설계 결과가 없습니다.');
 select(panels[0]?.modality);
 // 03: comparable route-level observations, not cross-modality molecular scores.
 const wrap=el(box($('pr-comparison'),'MODALITY COMPARISON'),'div',null,'pr-table'),table=el(wrap,'table');
 const head=el(el(table,'thead'),'tr'),tbody=el(table,'tbody');
 ['모달리티','실행','Validation','Critic','실행 방식','실험 검증','미평가 단계'].forEach(t=>el(head,'th',t));
 for(const p of panels){const row=el(tbody,'tr');[p.title,p.execution_status,p.validation,p.critic,mode(p),p.experimental_validation,(p.missing_steps||[]).join(' · ')||'미기록'].forEach(v=>el(row,'td',value(v)));}
 // 05: preserve recorded constraints, label suggested wet-lab work explicitly.
 const wet={SMALL_MOLECULE:['생화학적 활성·결합 평가','세포 내 target engagement 및 선택성 평가'],DE_NOVO_BINDER:['발현·정제 및 응집 평가','SPR/BLI 결합 및 비표적 결합 평가'],ADC:['내재화 및 정상조직 counter-screen','Linker·payload·DAR 및 접합체 세포독성 평가'],DEGRADER:['DC50/Dmax 및 E3·proteasome 의존성 확인','Proteomics 선택성 평가']};
 for(const p of panels){
  const b=box($('pr-validation'),p.title),content=el(b,'div',null,'qbody');
  note(content,'실험 검증: '+value(p.experimental_validation));
  if(p.blockers?.length){el(content,'b','판정 근거 · 제한 사항');const ul=el(content,'ul');p.blockers.forEach(t=>el(ul,'li',t));}
  if(p.missing_steps?.length)note(content,'미평가: '+p.missing_steps.join(' · '));
  if(Object.keys(p.computational_validation||{}).length)raw(content,'계산 검증 규칙 · 임계값',p.computational_validation);
  el(content,'b','후속 실험 제안 · 미실행');const ol=el(content,'ol');(wet[p.modality]||[]).forEach(t=>el(ol,'li',t));
 }
 // 06: immutable source records, including reused-run provenance.
 note($('pr-files'),'실행 ID · '+d.run_id);
 for(const p of panels)raw($('pr-files'),p.title+' · '+p.run_id,{execution_mode:p.execution_mode,stages:p.stages,events:p.events});
 raw($('pr-files'),'전체 실행 기록 · 판정 원본',d);
 $('pr-print').onclick=()=>window.print();
 if(d.running&&!window.__COMPETITION_RUN__){clearTimeout(window.omicraftReportTimer);window.omicraftReportTimer=setTimeout(()=>location.reload(),10000);}
}
