/*
 * The Exomiser - A tool to annotate and prioritize genomic variants
 *
 * Copyright (c) 2016-2020 Queen Mary University of London.
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

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package org.monarchinitiative.exomiser.core.writers;

import de.charite.compbio.jannovar.mendel.ModeOfInheritance;
import org.monarchinitiative.exomiser.core.analysis.AnalysisResults;

/**
 * 
 * @author Jules Jacobsen <jules.jacobsen@sanger.ac.uk>
 */
public interface ResultsWriter {
    
    /**
     * Writes the result data out to the file specified in the OutputSettings object for the specified mode of inheritance.
     *
     * @param modeOfInheritance
     * @param analysisResults
     * @param settings
     * @deprecated Use the {@link ResultsWriter#writeFile(AnalysisResults, OutputSettings)} method instead.
     */
    @Deprecated(forRemoval = true)
    default void writeFile(ModeOfInheritance modeOfInheritance, AnalysisResults analysisResults, OutputSettings settings) {
        writeFile(analysisResults, settings);
    }

    /**
     * Writes the result data out to the file specified in the OutputSettings object.
     *
     * @param analysisResults
     * @param outputSettings
     * @since 13.2.0
     */
    void writeFile(AnalysisResults analysisResults, OutputSettings outputSettings);

    /**
     * Writes the result data out to a String for the specified mode of inheritance.
     *
     * @param modeOfInheritance
     * @param analysisResults
     * @param settings
     * @return The string value of the analysis results for the given inheritance model
     * @deprecated use the {@link ResultsWriter#writeString(AnalysisResults, OutputSettings)} method instead.
     */
    @Deprecated(forRemoval = true)
    default String writeString(ModeOfInheritance modeOfInheritance, AnalysisResults analysisResults, OutputSettings settings) {
        return writeString(analysisResults, settings);
    }

    /**
     * Writes the result data out to a String, using the criteria provided by the {@link OutputSettings}.
     *
     * @param analysisResults
     * @param outputSettings
     * @return The string value of the analysis results for the given outputSettings.
     * @since 13.2.0
     */
    String writeString(AnalysisResults analysisResults, OutputSettings outputSettings);

}
