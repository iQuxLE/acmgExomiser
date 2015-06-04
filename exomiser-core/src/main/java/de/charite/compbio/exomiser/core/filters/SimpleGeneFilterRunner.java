/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package de.charite.compbio.exomiser.core.filters;

import de.charite.compbio.exomiser.core.model.Filterable;
import de.charite.compbio.exomiser.core.model.Gene;
import de.charite.compbio.exomiser.core.model.VariantEvaluation;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 *
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public class SimpleGeneFilterRunner implements FilterRunner<GeneFilter, Gene> {

    private static final Logger logger = LoggerFactory.getLogger(SimpleGeneFilterRunner.class);

    @Override
    public List<Gene> run(List<GeneFilter> filters, List<Gene> genes) {
        logger.info("Filtering {} genes using non-destructive simple filtering", genes.size());
        for (Gene gene : genes) {
            //Gene filtering needs to happen after variant filtering and only on genes which have passed the variant filtering steps
            //TODO: does this really have to be the case???
            if (gene.passedFilters()) {
                runAllFiltersOverGene(filters, gene);
            }
        }
        logger.info("Ran {} filters over {} genes using non-destructive simple filtering.", getFilterTypes(filters), genes.size());
        return genes;
    }

    @Override
    public List<Gene> run(GeneFilter filter, List<Gene> genes) {
        for (Gene gene : genes) {
            if (gene.passedFilters()) {
                runFilterAndAddResult(filter, gene);
            }
        }
        return genes;
    }

    private void runAllFiltersOverGene(List<GeneFilter> filters, Gene gene) {
        for (Filter filter : filters) {
            FilterResult filterResult = runFilterAndAddResult(filter, gene);
            //TODO: should this always be the case?
            addFilterResultToVariants(filterResult, gene.getVariantEvaluations());
        }
    }

    private void addFilterResultToVariants(FilterResult filterResult, List<VariantEvaluation> variantEvaluations) {
        for (VariantEvaluation variantEvaluation : variantEvaluations) {
            variantEvaluation.addFilterResult(filterResult);
        }
    }

    private FilterResult runFilterAndAddResult(Filter filter, Filterable filterable) {
        FilterResult filterResult = filter.runFilter(filterable);
        if (filterResult.getResultStatus() != FilterResultStatus.NOT_RUN) {
            filterable.addFilterResult(filterResult);
        }
        return filterResult;
    }

    private Set<FilterType> getFilterTypes(List<GeneFilter> filters) {
        Set<FilterType> filtersRun = new LinkedHashSet<>();
        for (Filter filter : filters) {
            filtersRun.add(filter.getFilterType());
        }
        return filtersRun;
    }

}
