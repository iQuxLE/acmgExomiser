/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2021 Queen Mary University of London.
 * Copyright (c) 2012-2016 Charité Universitätsmedizin Berlin and Genome Research Ltd.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.monarchinitiative.exomiser.autoconfigure.genome;

import de.charite.compbio.jannovar.data.JannovarData;
import org.h2.mvstore.MVStore;
import org.monarchinitiative.exomiser.core.genome.GenomeAnalysisService;
import org.monarchinitiative.exomiser.core.genome.GenomeDataService;
import org.monarchinitiative.exomiser.core.genome.VariantAnnotator;
import org.monarchinitiative.exomiser.core.genome.VariantDataService;
import org.monarchinitiative.exomiser.core.genome.dao.AllelePropertiesDao;
import org.monarchinitiative.exomiser.core.genome.dao.FrequencyDao;
import org.monarchinitiative.exomiser.core.genome.dao.PathogenicityDao;
import org.monarchinitiative.exomiser.core.genome.dao.VariantWhiteList;

/**
 * Interface to mark classes from the {@link org.monarchinitiative.exomiser.core.genome} package which need to have
 * unique beans exposed to Spring in order that they are intercepted for caching.
 *
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public interface GenomeAnalysisServiceConfiguration {

    JannovarData jannovarData();

    MVStore mvStore();

    VariantAnnotator variantAnnotator();

    GenomeDataService genomeDataService();

    VariantDataService variantDataService();

    VariantWhiteList variantWhiteList();

    //The classes below require Spring to managed the caching
    GenomeAnalysisService genomeAnalysisService();

    AllelePropertiesDao allelePropertiesDao();

    FrequencyDao localFrequencyDao();

    PathogenicityDao remmDao();

    PathogenicityDao caddDao();

    PathogenicityDao testPathScoreDao();

    // Structural Variant Dao
    FrequencyDao svFrequencyDao();

    PathogenicityDao svPathogenicityDao();
}
