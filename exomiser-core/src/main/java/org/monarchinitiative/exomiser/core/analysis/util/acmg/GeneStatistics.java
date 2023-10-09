package org.monarchinitiative.exomiser.core.analysis.util.acmg;

import org.monarchinitiative.exomiser.core.genome.dao.ClinVarDao;
import org.monarchinitiative.exomiser.core.model.VariantEvaluation;
import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;

import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

/*
do we have Gene range query? for Intervl of gene?
--> ClinVarDao Interval function to retrieve the data.
- distanceToNearestGene ?
- geneConstraints?
-ClinVarParser does not parse into mvStore? Where is that? also add GeneSymbol and then build Map in MvStore for each Gene
--> GeneMap (see below

- Gene holds list of VariantEvaluations that affect this Gene)

- make GeneEvaluation (object instead VariantEvalaution)

- hmm easiest w
 */

/*
    data class for geneStatistics
    = give geneName when open
    = get example:
    stop_gained {P=1,       frameshift: { P=1,
                 LP=2,                    LP=2,
                 B=3,                     B=3,
                 LB=4,                    LB=4,
                 ...}                       ....}
 */

public class GeneStatistics {

    private String geneName;

    private VariantEvaluation variantEvaluation;

    private ClinVarData clinVarData;
    private ClinVarDao clinVarDao;
    private Map<String, Map<ClinVarData.ClinSig, Integer>> consequenceClinicalSignificanceCounts;

    private GeneStatistics(String geneName, Map<String, Map<ClinVarData.ClinSig, Integer>> consequenceClinicalSignificanceCounts, ClinVarData clinVarData, ClinVarDao clinVarDao){
        this.geneName = geneName;
        this.consequenceClinicalSignificanceCounts = consequenceClinicalSignificanceCounts;
        this.clinVarData = Objects.requireNonNull(clinVarData);
        this.clinVarDao = Objects.requireNonNull(clinVarDao);

    }

    public Map<String, Map<ClinVarData.ClinSig, Integer>> getConsequenceClinicalSignificanceCounts(){
        return consequenceClinicalSignificanceCounts;
    }

    public String getGeneName(){
        return geneName;
    }

    public void countConsequence(String consequence, ClinVarData.ClinSig clinicalSignificance) {
     consequenceClinicalSignificanceCounts.putIfAbsent(consequence, new HashMap<>());
     Map<ClinVarData.ClinSig, Integer> clinicalSignificanceCounts = consequenceClinicalSignificanceCounts.get(consequence);

    int count = clinicalSignificanceCounts.getOrDefault(clinicalSignificance, 0);
    clinicalSignificanceCounts.put(clinicalSignificance, count+1);

    //countConsequence(stop_gained

    }

}
