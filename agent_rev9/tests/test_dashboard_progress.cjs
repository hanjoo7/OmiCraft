// Replay production dashboard events without network, inference, or a browser.
const fs=require('node:fs'),path=require('node:path'),vm=require('node:vm'),assert=require('node:assert/strict');
const source=fs.readFileSync(path.join(__dirname,'../assets/dashboard.js'),'utf8');
const dom=new Map();
const context=vm.createContext({document:{getElementById(id){if(!dom.has(id))dom.set(id,{});return dom.get(id);}},console});
vm.runInContext(source.slice(0,source.indexOf("$('historySelect').onchange=")),context);
vm.runInContext("appendLog=()=>{};reportLink=()=>{};loadHistory=()=>{};loadChoices=()=>{};loadStructures=()=>{};",context);
function run(code){return vm.runInContext(code,context);}
function reset(){run("resetDetailedProgress();state.steps=emptySteps();state.report='upstream_test';state.running=true;");}
function event(e){run(`consume(${JSON.stringify(e)})`);}
function status(name){return run(`detailState[${JSON.stringify(name)}].status`);}
reset();
assert.equal(status('assess_modalities'),'waiting');
event({type:'node',node:'user_gate',gate_status:'ALLOWED',messages:[],errors:[]});
assert.equal(status('assess_modalities'),'done');assert.equal(status('user_selection_gate_node'),'done');
event({type:'node',node:'design',execution_status:'COMPLETED'});
assert.equal(status('automatic_design_node'),'done');assert.equal(run('state.steps.design.status'),'done');
for(const [gate,expected] of [['HOLD','blocked'],['PENDING','pending'],['BLOCKED_REJECT','blocked'],[null,'pending']]){
 reset();event({type:'node',node:'user_gate',gate_status:gate});
 assert.equal(status('assess_modalities'),'done');assert.equal(status('user_selection_gate_node'),expected);
 event({type:'done',execution_status:'COMPLETED'});
 assert.equal(status('user_selection_gate_node'),expected);assert.equal(status('automatic_design_node'),'skipped');
}
for(const [execution,expected] of [['PARTIAL','partial'],['BLOCKED','blocked'],['FAILED','failed'],['NOT_RUN','skipped'],['CANCELLED','cancelled']]){
 reset();event({type:'node',node:'design',execution_status:execution});event({type:'done',execution_status:'COMPLETED'});
 assert.equal(status('automatic_design_node'),expected);assert.equal(run('state.steps.design.status'),expected);
}
reset();event({type:'progress',node:'therapeutic_design'});
event({type:'node',node:'adc/prediction',execution_status:'COMPLETED'});
assert.equal(status('automatic_design_node'),'running');
event({type:'node',node:'automatic_design_node',execution_status:'PARTIAL'});
assert.equal(status('automatic_design_node'),'partial');
console.log('Dashboard replay passed: completed, held, pending, skipped, partial, cancelled, failed, and child-stage events.');
const fixture={run_id:'upstream_latest',summary:{candidates:[{gene:'MMP7',decision:'ADVANCE'},{gene:'COL9A3',decision:'ADVANCE'},{gene:'HIDDEN',decision:'ADVANCE'}]},candidate_display:{hidden_genes:['HIDDEN']},design_report:{panels:[{gene:'MMP7',execution_status:'COMPLETED',validation:'ADVANCE',critic:'ADVANCE',is_mock:false},{gene:'COL9A3',execution_status:'COMPLETED',validation:'REJECT',critic:'REJECT'}]}};
context.targetFixture=fixture;
assert.equal(run('targetsFromReport(targetFixture).filter(t=>t.designPassed).map(t=>t.gene).join()'),'MMP7');
assert.equal(run('targetsFromReport(targetFixture).length'),2);
assert.equal(run("passedDesign({execution_status:'COMPLETED',validation:'ADVANCE',critic:'HOLD'})"),false);
assert.equal(run("passedDesign({execution_status:'COMPLETED',validation:'ADVANCE',critic:'ADVANCE',is_mock:true})"),false);
console.log('Latest target result checks passed: qualification and design decisions remain distinct.');
