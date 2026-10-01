function showStructureFallback(models, reason) {
  const host = document.getElementById('viewer');
  host.replaceChildren();
  const canvas = document.createElement('canvas');
  canvas.style.cssText="display:block;width:100%;height:100%;touch-action:none";
  host.style.overflow="hidden";
  host.append(canvas);
  const ctx = canvas.getContext('2d');
  if (!ctx) throw new Error('Canvas unavailable');
  const colors = ['#68b9ff', '#edb55c', '#b193ef', '#66d6a3'];
  const groups = [];
  for (const [text, format] of models) {
    if (format === 'pdb') {
      const chains = new Map(), ligands = new Map();
      for (const line of text.split('\n')) {
        const atom=line.slice(12,16).trim(), residue=line.slice(17,20).trim();
        const protein=line.startsWith('ATOM') && atom==='CA';
        const ligand=line.startsWith('HETATM') && !['HOH','WAT','DOD'].includes(residue);
        if(!protein && !ligand) continue;
        const point=[Number(line.slice(30,38)),Number(line.slice(38,46)),Number(line.slice(46,54))];
        if(!point.every(Number.isFinite)) continue;
        const collection=protein?chains:ligands;
        const key=protein?line.slice(21,22):line.slice(17,27);
        if(!collection.has(key)) collection.set(key,[]);
        collection.get(key).push(point);
      }
      for (const points of chains.values()) groups.push({points, trace:true});
      for (const points of ligands.values()) groups.push({points, trace:false});
    } else if (format === 'sdf') {
      const lines = text.split('\n');
      const count = Number((lines[3] || '').slice(0,3));
      const points = lines.slice(4,4+count).map(line=>[Number(line.slice(0,10)),Number(line.slice(10,20)),Number(line.slice(20,30))]);
      groups.push({points, trace:false});
    }
  }
  const points = groups.flatMap(g=>g.points).filter(p=>p.every(Number.isFinite));
  if (!points.length) throw new Error('No structure coordinates');
  const center = [0,1,2].map(i=>points.reduce((s,p)=>s+p[i],0)/points.length);
  const radius = Math.max(...points.map(p=>Math.hypot(...p.map((x,i)=>x-center[i]))), 1);
  let yaw=0, pitch=0, zoom=1, drag=null;
  function draw() {
    canvas.width=host.clientWidth; canvas.height=host.clientHeight;
    ctx.fillStyle='#0d1117';ctx.fillRect(0,0,canvas.width,canvas.height);
    const scale=Math.min(canvas.width,canvas.height)*.42/radius*zoom;
    function project(p){
      const [x,y,z]=p.map((v,i)=>v-center[i]);
      const a=x*Math.cos(yaw)+z*Math.sin(yaw),b=z*Math.cos(yaw)-x*Math.sin(yaw);
      return [canvas.width/2+a*scale,canvas.height/2+(y*Math.cos(pitch)-b*Math.sin(pitch))*scale];
    }
    groups.forEach((g,i)=>{
      ctx.strokeStyle=ctx.fillStyle=colors[i%colors.length];ctx.lineWidth=g.trace?1.8:3;
      ctx.beginPath();g.points.forEach((p,j)=>{const [x,y]=project(p);if(g.trace){if(j)ctx.lineTo(x,y);else ctx.moveTo(x,y);}else{ctx.moveTo(x+3,y);ctx.arc(x,y,3,0,Math.PI*2);}});
      if(g.trace)ctx.stroke();else ctx.fill();
    });
    ctx.fillStyle='#c1cad4';ctx.font='12px system-ui';
    ctx.fillText('Coordinate view · CA traces / ligand atoms · drag to rotate, scroll to zoom',12,24);
  }
  canvas.onpointerdown=e=>{drag=[e.clientX,e.clientY];canvas.setPointerCapture(e.pointerId);};
  canvas.onpointermove=e=>{if(drag){yaw+=(e.clientX-drag[0])*.01;pitch+=(e.clientY-drag[1])*.01;drag=[e.clientX,e.clientY];draw();}};
  canvas.onpointerup=canvas.onpointercancel=canvas.onlostpointercapture=()=>{drag=null;};
  canvas.onwheel=e=>{e.preventDefault();zoom=Math.max(.2,Math.min(8,zoom*Math.exp(-e.deltaY*.001)));draw();};
  new ResizeObserver(draw).observe(host);draw();
  const notice=document.createElement('small');notice.textContent=' · WebGL unavailable; showing coordinates (no molecular surface).';
  document.querySelector('header').append(notice);
  window.viewerMode='canvas_coordinates';window.viewerWarning=String(reason);window.viewerError=null;window.viewerReady=true;
}
