function renderPipelineDesign(container, report, options={}){
 const el=(parent,tag,text,cls)=>{const n=document.createElement(tag);if(text!=null)n.textContent=text;if(cls)n.className=cls;parent.append(n);return n;};
 const host=el(container,'section',null,'integrated-report pipeline-design');host.id='pipeline-design-results';
 if(!options.referenceOnly){
 el(host,'h2',options.gene?options.gene+' · 설계 결과':'Design Results');
 el(host,'p',`Workflow ${report.workflow_status} · Report ${report.report_status}`,'ir-muted');
 const metrics=el(host,'div',null,'ir-metrics');
 [['Targets evaluated',report.evaluated_targets],['Candidates generated',report.candidate_count],['Routes completed',report.completed_routes],['Routes passed',report.passing_routes]].forEach(([label,value])=>{
  const box=el(metrics,'div');el(box,'small',label);el(box,'strong',value);
 });
 if(report.candidate_search?.requested_limit)el(host,'p',`Candidate search limit · ${report.candidate_search.requested_limit} targets`,'ir-muted');
 el(host,'p',options.gene?'선택한 표적의 원래 실행 결과입니다. 후속 검증 결과는 아래에 별도로 표시합니다.':'Current run · 아래 수치와 판정은 이번 실행의 Design 결과입니다. 표적 수는 생성된 분자 후보 수와 다릅니다.','ir-muted');
 if(!report.panels.length)el(host,'p','No Design route was executed. Qualification results and recorded reasons remain available above.');
 else if(!report.passing_routes)el(host,'p','No candidate passed in this run. Completed computations, rejected candidates and blocked routes are retained below.');
 const preferred=report.panels.find(designRoutePassed)||report.panels.find(p=>p.execution_status==='COMPLETED');
 for(const p of report.panels){
  const detail=el(host,'details',null,'pipeline-design-route');detail.open=designRoutePassed(p)||p===preferred;if(designRoutePassed(p))detail.classList.add('ir-passed-route');
  el(detail,'summary',`${p.gene} · ${p.title} · ${p.execution_status} / ${p.validation}${designRoutePassed(p)?' · ★ 계산 검증 통과':''}`);
  const section=el(detail,'div',null,'ir-section');renderDesignPanel(section,{...p,title:p.gene+' · '+p.title},options);
 }
 if(report.followups?.length){
  const followups=el(host,'section',null,'pipeline-reference');
  el(followups,'h2','후속 설계 · 검증 결과');
  el(followups,'p','동일한 분석 실행에서 선택된 표적의 후속 실험입니다. 원래 실행의 통과 수에 합산하지 않습니다. 기존 후보 재사용은 실행 방식과 출처를 함께 표시합니다.','ir-muted');
  const preferred=report.followups.find(p=>designRoutePassed(p)&&!p.reuse);
  for(const p of report.followups){
   const detail=el(followups,'details');detail.open=designRoutePassed(p)||p===preferred;if(designRoutePassed(p))detail.classList.add('ir-passed-route');
   el(detail,'summary',`${p.gene} · ${p.title} · ${p.execution_status} / ${p.validation}${designRoutePassed(p)?' · ★ 계산 검증 통과':''} · ${p.execution_mode}`);
   if(p.reuse){el(detail,'p','기존 결과 재사용 · '+p.reuse.source_run_id,'ir-muted');}
   renderDesignPanel(el(detail,'div',null,'ir-section'),p,options);
  }
 }
 }
 if(report.reference){
  const reference=el(host,'section',null,'pipeline-reference');
  el(reference,'h2','Previously recorded · ERBB2');
  el(reference,'p','이전 실행의 실제 계산 결과입니다. 이번 실행에서 새로 생성된 후보나 통과 수에 포함하지 않습니다.','ir-muted');
  if(options.downloads!==false){const link=el(reference,'a','Open ERBB2 integrated report ↗');link.href=report.reference.report_url;}
  el(reference,'p',report.reference.run_id,'ir-muted');
  const preferred=report.reference.panels.find(designRoutePassed);
  for(const p of report.reference.panels){
   const detail=el(reference,'details');detail.open=designRoutePassed(p)||p===preferred;if(designRoutePassed(p))detail.classList.add('ir-passed-route');
   el(detail,'summary',`${p.title} · ${p.execution_status} / ${p.validation}${designRoutePassed(p)?' · ★ 계산 검증 통과':''} · Previous run`);
   renderDesignPanel(el(detail,'div',null,'ir-section'),p,options);
  }
 }
 if(report.reference_warning)el(host,'p','Previous report unavailable: '+report.reference_warning,'ir-muted');
 return host;
}
