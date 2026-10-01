const targetStructureCatalog={items:[],status:'loading',pending:null};
function renderTargetStructures(row,gene){
 row.dataset.gene=String(gene||'').toUpperCase();
 row.querySelector('.structure-links')?.remove();
 const links=document.createElement('div');links.className='structure-links';row.append(links);
 const matches=targetStructureCatalog.items.filter(asset=>String(asset.gene).toUpperCase()===row.dataset.gene);
 if(!matches.length){
  const button=document.createElement('button');button.type='button';button.disabled=true;button.className='struct-link';
  const status=targetStructureCatalog.status;
  button.textContent='🔬 Structure · '+(status==='loading'?'Loading':status==='error'?'Unavailable':'Not available');
  button.title=status==='ready'?'No registered structure for '+row.dataset.gene:'Structure catalog '+status;
  links.append(button);return;
 }
 for(const asset of matches){
  const source=String(asset.source||'');let label;
  if(asset.source_type==='alphafold_db'||/AlphaFold DB/i.test(source))label='AlphaFold DB Structure';
  else if(/\b(?:AF3|AlphaFold[ -]?3)\b/i.test(source))label='AF3 Structure · '+asset.label;
  else if(asset.source_type==='pdb'||/^PDB /i.test(source))label='PDB Structure · '+(asset.pdb_id||asset.label);
  else label='Structure · '+asset.label;
  const link=document.createElement('a');link.className='struct-link';link.textContent='🔬 '+label;
  link.href='/api/structures/'+encodeURIComponent(asset.id)+'/viewer';link.target='_blank';link.rel='noopener';
  link.title=asset.label+' — '+source;links.append(link);
 }
}
function refreshTargetStructures(){
 if(targetStructureCatalog.pending)return targetStructureCatalog.pending;
 targetStructureCatalog.pending=(async()=>{
  try{
   const response=await fetch('/api/structures',{cache:'no-store'});
   if(!response.ok)throw new Error('Structure catalog HTTP '+response.status);
   targetStructureCatalog.items=(await response.json()).filter(asset=>asset.kind==='protein');
   targetStructureCatalog.status='ready';
  }catch(error){targetStructureCatalog.status='error';console.warn('Structure catalog:',error.message);}
  finally{
   targetStructureCatalog.pending=null;
   document.querySelectorAll('#tgtL .tgt[data-gene]').forEach(row=>renderTargetStructures(row,row.dataset.gene));
  }
 })();return targetStructureCatalog.pending;
}
