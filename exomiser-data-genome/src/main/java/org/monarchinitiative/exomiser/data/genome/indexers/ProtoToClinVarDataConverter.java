package org.monarchinitiative.exomiser.data.genome.indexers;

import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
import org.monarchinitiative.exomiser.core.proto.AlleleProto;

import java.util.HashSet;
import java.util.Set;

public class ProtoToClinVarDataConverter {

    /**
     * Converts a protobuf ClinVar object into a ClinVarData object.
     */
    public static ClinVarData fromProtoClinVar(AlleleProto.ClinVar clinVar) {
        ClinVarData.Builder builder = ClinVarData.builder();

        builder.alleleId(clinVar.getAlleleId())
                .primaryInterpretation(fromProtoClinSig(clinVar.getPrimaryInterpretation()));

        // If your ClinVarData supports setting secondary interpretations, you can do something like this:
        Set<ClinVarData.ClinSig> secondaryInterpretationsSet = new HashSet<>();

        for (AlleleProto.ClinVar.ClinSig clinSig : clinVar.getSecondaryInterpretationsList()) {
            secondaryInterpretationsSet.add(fromProtoClinSig(clinSig));
        }
        builder.secondaryInterpretations(secondaryInterpretationsSet);

        // Add other fields as required, e.g., reviewStatus, includedAlleles

        return builder.build();
    }

    /**
     * Converts a protobuf ClinSig enum value into a ClinVarData.ClinSig enum value.
     */
    public static ClinVarData.ClinSig fromProtoClinSig(AlleleProto.ClinVar.ClinSig clinSig) {
        switch (clinSig) {
            case BENIGN:
                return ClinVarData.ClinSig.BENIGN;
            case BENIGN_OR_LIKELY_BENIGN:
                return ClinVarData.ClinSig.BENIGN_OR_LIKELY_BENIGN;
            case LIKELY_BENIGN:
                return ClinVarData.ClinSig.LIKELY_BENIGN;
            case UNCERTAIN_SIGNIFICANCE:
                return ClinVarData.ClinSig.UNCERTAIN_SIGNIFICANCE;
            case LIKELY_PATHOGENIC:
                return ClinVarData.ClinSig.LIKELY_PATHOGENIC;
            case PATHOGENIC_OR_LIKELY_PATHOGENIC:
                return ClinVarData.ClinSig.PATHOGENIC_OR_LIKELY_PATHOGENIC;
            case PATHOGENIC:
                return ClinVarData.ClinSig.PATHOGENIC;
            default:
                throw new IllegalArgumentException(clinSig + " not a recognised value");
        }
    }
}

