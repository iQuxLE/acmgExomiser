/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.core;

import de.charite.compbio.exomiser.core.factories.VariantDataService;
import de.charite.compbio.exomiser.core.filters.FilterFactory;
import de.charite.compbio.exomiser.core.filters.GeneFilter;
import de.charite.compbio.exomiser.core.filters.VariantFilter;
import de.charite.compbio.exomiser.core.model.SampleData;
import de.charite.compbio.exomiser.core.prioritisers.Prioritiser;
import de.charite.compbio.exomiser.core.prioritisers.PriorityFactory;
import de.charite.compbio.exomiser.core.prioritisers.PriorityType;
import java.util.ArrayList;
import java.util.List;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Main class for analysing variant data. This will orchestrate the set-up of
 * Filters and Priotitisers according to the supplied settings and then apply
 * them to the data.
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class Exomiser {

    private static final Logger logger = LoggerFactory.getLogger(Exomiser.class);

    private final PriorityFactory prioritiserFactory;
    private final VariantDataService variantDataService;

    public Exomiser(VariantDataService variantDataService, PriorityFactory prioritiserFactory) {
        this.variantDataService = variantDataService;
        this.prioritiserFactory = prioritiserFactory;
    }

    //TODO: this could be made more programmatically configuarable by also allowing
    //the Filters and Prioritisers be passsed in / set so as to remove the dependency on ExomiserSettings
    //will happen, but not through here - Analysis, AnalysisRunner and some missing bits will do that.
    public Analysis analyse(SampleData sampleData, ExomiserSettings exomiserSettings) {
        
        Analysis analysis = setUpExomiserAnalysis(sampleData, exomiserSettings);
        //TODO: make run the functional entry point, remove all the other methods or push them into AnalysisRunner
        //AnalysisRunner is kind of the new Exomiser - exomiser is now just a brand and a set of immutable steps 
        //from setUpExomiserAnalysis played out by the AnalysisRunner and stored in an Analysis
        logger.info("STARTING ANALYSIS");
        AnalysisRunner analysisRunner = new AnalysisRunner(variantDataService, exomiserSettings);
        analysisRunner.runAnalysis(analysis);
        logger.info("FINISHED ANALYSIS");
        return analysis;
    }

    /**
     * Sets up an Analysis using the published Exomiser algorithm where variants
     * are filtered first, then genes and finally the prioritisers are run.
     *
     * @param sampleData
     * @param exomiserSettings
     * @return
     */
    private Analysis setUpExomiserAnalysis(SampleData sampleData, ExomiserSettings exomiserSettings) {
        logger.info("SETTING-UP ANALYSIS");
        Analysis analysis = new Analysis(sampleData, exomiserSettings);
        
        //don't change the order here - variants should ALWAYS be filtered before
        addVariantFilters(exomiserSettings, analysis);
        //genes have their inheritance modes filtered otherwise the inheritance mode will break leading to altered
        //predictions downstream as we only want to test the mode for candidate variants.
        addGeneFilters(exomiserSettings, analysis);
        //inheritance modes are needed even if we don't have an inheritance gene filter set as the OMIM prioritiser relies on it
        //Prioritisers should ALWAYS run last.
        addPrioritisers(exomiserSettings, analysis);

        return analysis;
    }

    private List<VariantFilter> addVariantFilters(ExomiserSettings exomiserSettings, Analysis analysis) {
        FilterFactory filterFactory = new FilterFactory();
        logger.info("MAKING VARIANT FILTERS");
        List<VariantFilter> variantFilters = filterFactory.makeVariantFilters(exomiserSettings);
        for (VariantFilter variantFilter : variantFilters) {
            analysis.addStep(variantFilter);
        }
        return variantFilters;
    }

    private List<GeneFilter> addGeneFilters(ExomiserSettings exomiserSettings, Analysis analysis) {
        FilterFactory filterFactory = new FilterFactory();
        logger.info("MAKING GENE FILTERS");
        List<GeneFilter> geneFilters = filterFactory.makeGeneFilters(exomiserSettings);
        for (GeneFilter geneFilter : geneFilters) {
            analysis.addStep(geneFilter);
        }
        return geneFilters;
    }

    private List<Prioritiser> addPrioritisers(ExomiserSettings exomiserSettings, Analysis analysis) {
        logger.info("MAKING PRIORITISERS");
        List<Prioritiser> prioritisers = makePrioritisers(exomiserSettings);
        for (Prioritiser prioritiser : prioritisers) {
            analysis.addStep(prioritiser);
        }
        return prioritisers;
    }

    private List<Prioritiser> makePrioritisers(ExomiserSettings exomiserSettings) {
        List<Prioritiser> prioritisers = new ArrayList<>();

        PriorityType prioritiserType = exomiserSettings.getPrioritiserType();
        if (prioritiserType == PriorityType.NONE || prioritiserType == PriorityType.NOT_SET) {
            return prioritisers;
        }

        //always run OMIM unless the user specified what they really don't want to run any prioritisers
        Prioritiser omimPrioritiser = prioritiserFactory.makePrioritiser(PriorityType.OMIM_PRIORITY, exomiserSettings);
        prioritisers.add(omimPrioritiser);
        //don't add OMIM prioritiser twice to the list
        if (prioritiserType != PriorityType.OMIM_PRIORITY) {
            Prioritiser prioritiser = prioritiserFactory.makePrioritiser(prioritiserType, exomiserSettings);
            prioritisers.add(prioritiser);
        }
        return prioritisers;
    }

}
