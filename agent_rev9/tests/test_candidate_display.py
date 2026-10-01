import copy
import json
from types import SimpleNamespace

import pytest

from .. import candidate_display as display, target_qualification as qual, web_execution as web
from ..tier_system import assign_bio_confidence
from ..pipeline_report import assemble_pipeline_report


@pytest.mark.parametrize('genetic,somatic,tier', [(None,None,'B3'), (0,0,'B3'), (.1,None,'B2'), (None,.1,'B2'), (.1,.2,'B2')])
def test_ot_independent_support(genetic, somatic, tier):
    b, reasons = assign_bio_confidence({'dual_contrast':'dual_contrast_supported',
        'ot_genetic_association': genetic, 'ot_somatic_mutation': somatic})
    assert b == tier
    assert sum(r.startswith('OT ') for r in reasons) == sum(v is not None and v > 0 for v in (genetic, somatic))


def test_ot_can_reach_b1_but_never_override_veto():
    evidence = {'dual_contrast':'partial_evidence','has_functional_evidence':True,
                'ot_genetic_association':.1,'ot_somatic_mutation':.2}
    assert assign_bio_confidence(evidence)[0] == 'B1'
    assert assign_bio_confidence({**evidence,'has_safety_veto':True})[0] == 'B0'
    assert assign_bio_confidence({**evidence,'dual_contrast':'discordant'})[0] == 'B0'


@pytest.mark.parametrize('status,expected', [('FOUND','B2'),('NOT_FOUND','B3'),('ERROR','B3')])
def test_qualification_only_scores_found_ot(monkeypatch, status, expected):
    monkeypatch.setattr(qual,'query_uniprot',lambda *a:{'status':'NOT_FOUND'})
    monkeypatch.setattr(qual,'query_hpa_normal',lambda *a:{'status':'NOT_FOUND'})
    monkeypatch.setattr(qual,'query_chembl',lambda *a:{'status':'NOT_FOUND'})
    monkeypatch.setattr(qual,'query_open_targets',lambda *a:{'status':status,'genetic_association_score':.1,'somatic_mutation_score':.2})
    result=qual.qualify_gene('TEST',{'gene_id':'ENSG00000000001','dual_contrast':'partial_evidence'}, {}, {}, modalities=['DE_NOVO_BINDER'])
    assert result['bio_confidence'] == expected


def test_filters_share_lab_report_rules_without_mutating_records(tmp_path,monkeypatch):
    dataset=tmp_path/'dataset';(dataset/'depmap').mkdir(parents=True)
    (dataset/'depmap/common_essentials_computed.csv').write_text('gene\nESSENTIAL\n')
    dc=dataset/'tcga_brca/processed/deg_results';dc.mkdir(parents=True)
    (dc/'dual_contrast_classification.csv').write_text('gene_name,c1_direction\nOLD_DOWN,Down\nKEEP,Up\n')
    monkeypatch.setenv('OMICRAFT_DATASET_ROOT',str(dataset))
    root=tmp_path/'runs';folder=root/'upstream_test';plots=folder/'analysis/04_de_plots';plots.mkdir(parents=True)
    volcano=plots/'volcano_data.tsv'
    volcano.write_text('gene_symbol\tlog2FoldChange\tminus_log10_padj\tdirection\nNEW_DOWN\t-2\t4\tDown_in_TNBC\nKEEP\t2\t5\tUp_in_TNBC\n')
    genes=['ESSENTIAL','OLD_DOWN','NEW_DOWN','CDKN2A','KEEP']
    state={'r_analysis_result':{'output_dir':str(plots.parent)},'qualification_results':{g:{'tier_results':[{'modality':'DE_NOVO_BINDER','Tier':'2B'}]} for g in genes}}
    (folder/'pipeline_state.json').write_text(json.dumps(state))
    monkeypatch.setattr(web,'get_config',lambda:SimpleNamespace(data=SimpleNamespace(base_dir=str(root))))
    rows=web.design_candidates()
    assert [r['gene'] for r in rows if r['origin']=='upstream_analysis'] == ['KEEP']
    assert len([r for r in rows if r['origin']=='prior_reference_inputs']) == 4
    data={'run_id':folder.name,'summary':{'candidates':[{'gene':g} for g in genes]},'artifacts':[{'path':str(volcano)}]}
    original=copy.deepcopy(data)
    result=assemble_pipeline_report(data,root)
    assert set(result['candidate_display']['hidden_genes']) == set(genes)-{'KEEP'}
    assert result['summary']['candidates'] == original['summary']['candidates']
    assert data == original
    assert json.loads((folder/'pipeline_state.json').read_text()) == state
    assert display.visible_gene(' cdkn2a ') is False


def test_missing_filter_files_do_not_invent_exclusions(tmp_path,monkeypatch):
    monkeypatch.setenv('OMICRAFT_DATASET_ROOT',str(tmp_path))
    assert display.hidden_genes() == {'CDKN2A'}
