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

package org.monarchinitiative.exomiser.data.genome.model;


import org.monarchinitiative.exomiser.data.genome.indexers.Indexer;
import org.monarchinitiative.exomiser.data.genome.indexers.IndexingException;

import java.io.IOException;


/**
 * @author Jules Jacobsen <j.jacobsen@qmul.ac.uk>
 */
public interface OutputFileIndexingResource<T extends OutputLine> extends Resource<T> {

    Indexer<T> indexer();

    default void indexResource() {
        try (Indexer<T> indexer = indexer()) {
            indexer.index(this);
        } catch (IOException exception) {
            throw new IndexingException(exception);
        }
    }

}
