(() => {
  function simplifyDiagnostic() {
    if (document.body.id !== 'visual-diagnostic') return;
    ['issues', 'provenance'].forEach(id => {
      document.getElementById('tab-' + id)?.remove();
      document.getElementById('panel-' + id)?.remove();
    });
    document.getElementById('show-all')?.remove();
    document.querySelectorAll('#panel-small-molecule figure').forEach(figure => {
      if (/admet summary|endpoints/i.test(figure.querySelector('img')?.alt || '')) figure.remove();
    });
    const tabs = [...document.querySelectorAll('#visual-diagnostic [role=tab]')];
    const cards = [...document.querySelectorAll('#visual-diagnostic .vd-card')];
    const select = tab => {
      tabs.forEach(button => {
        const selected = button === tab;
        button.setAttribute('aria-selected', String(selected));
        document.getElementById(button.getAttribute('aria-controls')).hidden = !selected;
      });
      cards.forEach(card => card.setAttribute('aria-pressed', String(card.querySelector('h3')?.textContent === tab.textContent)));
    };
    tabs.forEach(tab => { tab.onclick = () => select(tab); });
    cards.forEach(card => {
      const tab = tabs.find(tab => tab.textContent === card.querySelector('h3')?.textContent);
      if (!tab) return;
      card.setAttribute('role', 'button');card.tabIndex = 0;
      card.setAttribute('aria-controls', tab.getAttribute('aria-controls'));
      card.onclick = () => { select(tab); tab.parentElement.scrollIntoView({block: 'start'}); };
      card.onkeydown = event => { if (event.key === 'Enter' || event.key === ' ') { event.preventDefault();card.click(); } };
    });
    const selected = tabs.find(tab => tab.getAttribute('aria-selected') === 'true') || tabs[0];
    if (selected) select(selected);
  }
  function styleContent() {
    simplifyDiagnostic();
    const host = document.querySelector('.report-body') || document.querySelector('.content');
    if (host && Array.from(host.children).some(node => /^H[123]$/.test(node.tagName))) {
      let card;
      for (const node of Array.from(host.children)) {
        if (node.tagName === 'STYLE' || node.tagName === 'SCRIPT') continue;
        if (!card || /^H[123]$/.test(node.tagName)) {
          card = document.createElement('section');
          card.className = 'report-card';
          host.insertBefore(card, node);
        }
        card.append(node);
      }
    }
    document.querySelectorAll('.report-body img, .content img, #visual-diagnostic main img').forEach(image => {
      if (image.closest('a') || !image.getAttribute('src')) return;
      const link = document.createElement('a');
      link.className = 'report-figure-link';
      link.href = image.getAttribute('src');
      link.target = '_blank';link.rel = 'noopener';link.title = 'Open full-size image';
      image.before(link);link.append(image);
    });
    document.querySelectorAll('.report-body table, .content table, #visual-diagnostic main table').forEach(table => {
      if (table.parentElement.matches('.report-table-scroll,.table-wrap,.scroll')) return;
      const wrapper = document.createElement('div');
      wrapper.className = 'report-table-scroll';
      wrapper.tabIndex = 0;
      wrapper.setAttribute('role', 'region');
      wrapper.setAttribute('aria-label', 'Report table');
      table.before(wrapper);
      wrapper.append(table);
    });
  }
  function mountReportTheme() {
    if (!document.getElementById('omicraft-report-brand')) {
      const template = document.createElement('template');
      template.innerHTML = __REPORT_BRAND__;
      document.body.prepend(template.content);
      document.getElementById('reportPrint').onclick = () => window.print();
    }
    styleContent();
  }
  const render = window.renderTherapeuticRun;
  if (render) window.renderTherapeuticRun = function (data) {
    const result = render(data);
    mountReportTheme();
    return result;
  };
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', mountReportTheme, {once: true});
  else mountReportTheme();
})();
