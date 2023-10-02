//package org.monarchinitiative.exomiser.core.genome;
//
//import org.monarchinitiative.exomiser.core.model.*;
//import org.monarchinitiative.exomiser.core.model.frequency.FrequencyData;
//import org.monarchinitiative.exomiser.core.model.frequency.FrequencySource;
//import org.monarchinitiative.exomiser.core.model.pathogenicity.ClinVarData;
//import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicityData;
//import org.monarchinitiative.exomiser.core.model.pathogenicity.PathogenicitySource;
//import org.monarchinitiative.svart.GenomicInterval;
//import org.monarchinitiative.svart.GenomicVariant;
//
//import java.util.List;
//import java.util.Map;
//import java.util.Set;
//
//public class TestGenomeAnalysisService implements GenomeAnalysisService {
//
//    private final GenomeAssembly genomeAssembly;
//
//    private final GenomeDataService genomeDataService;
//    private final VariantDataService variantDataService;
//    private final VariantAnnotator variantAnnotator;
//
//    public TestGenomeAnalysisService(Builder builder) {
//        this.genomeAssembly = builder.genomeAssembly;
//        this.genomeDataService = builder.genomeDataService;
//        this.variantDataService = builder.variantDataService;
//        this.variantAnnotator = builder.variantAnnotator;
//    }
//
//    @Override
//    public GenomeAssembly getGenomeAssembly() {
//        return null;
//    }
//
//    @Override
//    public VariantAnnotator getVariantAnnotator() {
//        return null;
//    }
//
//
//}
