/*

timbreID - A generic classification external.

Copyright 2009 William Brent

This file is part of timbreID.

timbreID is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

timbreID is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

** 0.7.0 fixed concatID method to give it an option to do wraparound or not. Updated ID method to report a desired number of matches (based on KNN setting), not just one. This will help with doing timbre blending in timbre space example, where we need the N closest matches. Changed write_text method so that if the database is normalized, the normalized data is printed out. Also works with attribute_range and attributeOrder so that you can print out a specific subset of attribute columns instead of the entire attribute list. Still need to update ID method to remove the dynamically assigned *votes memory.

** 0.6.0 changed all floats to t_floats, got rid of unused headers, removed underscores from variable names.  FIXED concatenative_id function's search center/neighborhood bug.  Added a similarity matrix function.  Added a "print" method to show internal settings.  Made even or odd neighborhood settings valid since it will sometimes be important to pick an exact range.  Odd settings are fine - with n=5, we search from the center -/+ 2.  With n=6, it be center -2/+3 because of a special case for searchStart.
** 0.5.1 was just to remove the "no such table" error from the NRT externs when it doesn't find the specified table at load time. also, the max window size for all the NRT externs was set to 131072 because this seems to be the maximum allowable size for mayer_realfft().
** 0.5 as part of the timbreID-0.5 update, this gets rid of unnecessary getbytes (0) calls in timbreID.c. does not yet use memset() or memcpy() where possible.  fixed lack of memory free for x->x_weights, x->x_attributeOrder, x->x_instanceFeatureLengths, x->x_normData in _free.
** fixed small memory mismanagement in all _make_filterbank functions for bark-based externs.
** [timbreID] is unchanged - only updating for addition of [featureAccum] and [binWrangler], which streamline time-evolving feature extraction
** [timbreID] is unchanged - only updating version for [cepstrum] and [cepstrum~] changes.
** part of the update acknowledging the change to one tID folder full of sources and binaries. This update also included filter power averaging in barkSpec, bfcc, and mfcc.
** removed cosine similarity distance metric
** 0.3.3 reflects the general package update, which includes specIrregularity (~), and fixes a bug in cepstrum (~). also, changed timbreID's default KNN setting to 1.  It should really only be anything different if clustering is done.
** 0.3.3 is based on the windows source, where the _read functions were changed to avoid using fdopen.  also corrected the use of a resizable array for semicolon in the _read_cluster_text function. also added static to all functions except _setup.
** 0.3.2 no changes to [timbreID]. Added RT and NRT spectral spread externs, plus mean and standard deviation externs for summarizing spectro-temporal features of various lengths
** 0.3.1G this gets around to adding the cluster read/write functions. F is cleaned up version of D. 0.3.1G does work with read/write clusters.
** 0.3.1F skips the different approach in E (didn't take that very far), and fixed D to work.  The problem in _compute_cluster was after clustering and in memory cleanup. Shifting elements of x->x_clusterMembers and x->x_clusterMemberLengths was done in the wrong order and memory was resized incorrectly. F cleans up D, and adds the cluster_write functions.
** 0.3.1D all that remains is cluster write/read.
** 0.3.1C extends 0.3.1B to shoot for making x->x_clusterMembers memory dynamically allocated as well.  C works up to and including unclustering - next will be the read/write functions.
** 0.3.1 finally makes feature database and cluster member memory dynamically allocated. Set up the _read and _write functions so that the header includes the entire instanceFeatureLengths array. Might as well keep all that information so that mixed instance databases won't be a problem in the future. BUT - write/read_text do not store/retrieve this information. More importantly, x->x_featureLength is simply set to the first element in x->x_instanceFeatureLengths in the _read function. So, the repercussions of mixed length features in the database has not been thought through. A common feature length is still assumed.
** 0.3.0 adds a method for finding the worst match, and two others for outputting the maxes and mins of all features in the database (as long as it has been normalized). Also - found a bug in the concat_id chain. timbreID_clear never set x->x_neighborhood back to 0, which caused crashes after training/clearing/training a few times.
** 0.2.9 removes some warnings about unitialized variables. the entire timbreID set of externs was updated on this date as well to include non real-time versions of all feature externs. 0.3.0 will finally implement dynamic memory allocation for feature database instances, rather than using MAXFEATURELENGTH for all (very wasteful).
** 0.2.8 feature externs now use mayer_realfft, switched to timbreID-help.pd rather than help-timbreID.pd for the help patch name.
** 0.2.7C Added clusterMembership function
** 0.2.7 Added ability to output a instance's feature list from the 4th outlet.
** 0.2.6 Added ability to write and read cluster instances in binary and text formats.
** 0.2.5 Added correlation as a distance metric. Fixed a bug in compute_cluster: it now throws an error if you ask for 0 clusters.  Also slightly tweaked the clustering algorithm itself.
** 0.2.4 Fixed major bug where timbreID_train wasn't updating the size of x->x_featureInput if the feature length changes. This caused crashes in timbreID_id after training with a feature length different than the default of 47 points. Fixed bug where requesting compute_cluster with a different number of clusters when already clustered causes crash. Fixed error in _write: j was used uninitialized. Distance metric choice is now used in compute_cluster. Previously, it was squared_euclid only. Added method for outputting cluster member lists, or the entire cluster list.
** 0.2.3 adds binary file output (.timid now default), and an option for text output.
** 0.2.2 adds MATLAB .mat file output.
** 0.2.0 adds cosine distance function (doesn't account for attribute normalization yet though!).
** 0.1.9 adds ARFF function to export data for WEKA.

** fixed bug where the distances of the final elements in a compute_order result were sometimes
larger than INT_MAX, so that they were all index 0 (init value). Using FLT_MAX instead now, which
should be big enough to be safe.

** Accounted for normalization flag in the compute_variance function (requiring a new timbreID_mean
function as well).

*/

#include "tIDLib.h"

static t_class* timbreID_class;

typedef enum
{
    euclidean = 0,
    taxi,
    correlation
} t_distMetric;

typedef struct _timbreID
{
    t_object x_obj;
    t_symbol* x_objSymbol;
    t_instance* x_instances;
    t_instance* x_classRefs;
    t_cluster* x_clusters;
    t_attributeData* x_attributeData;
    t_attributeIdx x_maxFeatureLength;
    t_attributeIdx x_minFeatureLength;
    t_instanceIdx x_numInstances;
    t_instanceIdx x_numClassRefs;
    t_instanceIdx x_numClusters;
    t_distMetric x_distMetric;
    t_instanceIdx x_k;
    t_bool x_outputKnnMatches;
    t_bool x_normalize;
    t_bool x_relativeOrdering;

    t_bool x_concatWrap;
    t_bool x_reorientFlag;
    t_instanceIdx x_neighborhood;
    t_instanceIdx x_searchCenter;
    t_instanceIdx x_prevMatch;
    t_instanceIdx x_maxMatches;
    t_bool x_stutterProtect;
    t_float x_jumpProb;

    t_attributeIdx x_attributeLo;
    t_attributeIdx x_attributeHi;
    t_canvas* x_canvas;
    t_outlet* x_id;
    t_outlet* x_nearestDist;
    t_outlet* x_confidence;
    t_outlet* x_listOut;
} t_timbreID;


/* -------- utility functions -------- */
static t_float timbreID_attributeMean (t_instanceIdx numRows, t_attributeIdx column, t_instance* instances, t_bool normalFlag, t_attributeData* attributeData)
{
    t_instanceIdx i;
    t_float avg, min, scalar;

    avg = 0.0;

    for (i = 0; i < numRows; i++)
    {
        // check that this attribute is in range for instance i
        if (column < instances[i].length)
        {
            if (normalFlag)
            {
                min = attributeData[column].normData.minVal;
                scalar = attributeData[column].normData.normScalar;

                avg += (instances[i].data[column] - min) * scalar;
            }
            else
                avg += instances[i].data[column];
        }
    }

    avg /= numRows;

    return (avg);
}


// this gives the distance between the current inputData and a given instance from x_instances. Same as timbreID_getDist() below, except more convenient for use in the _id() method because you don't have to fill 2 t_float buffers in order to use it.
static t_float timbreID_getInputDist (t_timbreID* x, t_instanceIdx instanceID)
{
    t_attributeIdx i, j, vecLen;
    t_float min, max, normScalar, dist;
    t_float* vec1Buffer;
    t_float* vec2Buffer;
    t_float* vecWeights;

    vecLen = x->x_attributeHi - x->x_attributeLo + 1;
    min = max = dist = 0.0;
    normScalar = 1.0;

    vec1Buffer = (t_float *)t_getbytes (vecLen * sizeof (t_float));
    vec2Buffer = (t_float *)t_getbytes (vecLen * sizeof (t_float));
    vecWeights = (t_float *)t_getbytes (vecLen * sizeof (t_float));

    // extract the right attributes and apply normalization if it's active
    for (i = x->x_attributeLo, j = 0; i <= x->x_attributeHi; i++, j++)
    {
        t_attributeIdx thisAttribute;
        thisAttribute = x->x_attributeData[i].order;

        if (thisAttribute > x->x_instances[instanceID].length)
        {
            pd_error (x, "%s: attribute %i out of range for database instance %i.", x->x_objSymbol->s_name, thisAttribute, instanceID);

            // free local memory before exiting
            t_freebytes (vec1Buffer, vecLen * sizeof (t_float));
            t_freebytes (vec2Buffer, vecLen * sizeof (t_float));
            t_freebytes (vecWeights, vecLen * sizeof (t_float));

            return (FLT_MAX);
        }

        vecWeights[j] = x->x_attributeData[thisAttribute].weight;

        if (x->x_normalize)
        {
            if (x->x_attributeData[thisAttribute].inputData < x->x_attributeData[thisAttribute].normData.minVal)
                min = x->x_attributeData[thisAttribute].inputData;
            else
                min = x->x_attributeData[thisAttribute].normData.minVal;

            if (x->x_attributeData[thisAttribute].inputData > x->x_attributeData[thisAttribute].normData.maxVal)
            {
                max = x->x_attributeData[thisAttribute].inputData;
                normScalar = 1.0 / (max - min);
            }
            else
            {
                max = x->x_attributeData[thisAttribute].normData.maxVal;
                normScalar = x->x_attributeData[thisAttribute].normData.normScalar;
            }

            vec1Buffer[j] = (x->x_attributeData[thisAttribute].inputData - min) * normScalar;
            vec2Buffer[j] = (x->x_instances[instanceID].data[thisAttribute] - min) * normScalar;
        }
        else
        {
            vec1Buffer[j] = x->x_attributeData[thisAttribute].inputData;
            vec2Buffer[j] = x->x_instances[instanceID].data[thisAttribute];
        }
    }

    switch (x->x_distMetric)
    {
        case euclidean:
            dist = tIDLib_euclidDist (vecLen, vec1Buffer, vec2Buffer, vecWeights, false);
            break;
        case taxi:
            dist = tIDLib_taxiDist (vecLen, vec1Buffer, vec2Buffer, vecWeights);
            break;
        case correlation:
            dist = tIDLib_corr (vecLen, vec1Buffer, vec2Buffer);
            // bash to the 0-2 range, then flip sign so that lower is better. this keeps things consistent with other distance metrics.
            dist += 1;
            dist *= -1;
            break;
        default:
            break;
    };

    // free local memory
    t_freebytes (vec1Buffer, vecLen * sizeof (t_float));
    t_freebytes (vec2Buffer, vecLen * sizeof (t_float));
    t_freebytes (vecWeights, vecLen * sizeof (t_float));

    return (dist);
}


// this gives the distance between two vectors stored in the .data part of t_instance arrays.
static t_float timbreID_getDist (t_timbreID* x, t_instance instance1, t_instance instance2)
{
    t_attributeIdx i, j, vecLen;
    t_float min, max, normScalar, dist;
    t_float* vec1Buffer;
    t_float* vec2Buffer;
    t_float* vecWeights;

    vecLen = x->x_attributeHi - x->x_attributeLo + 1;
    min = max = dist = 0.0;
    normScalar = 1.0;

    vec1Buffer = (t_float *)t_getbytes (vecLen * sizeof (t_float));
    vec2Buffer = (t_float *)t_getbytes (vecLen * sizeof (t_float));
    vecWeights = (t_float *)t_getbytes (vecLen * sizeof (t_float));

    // extract the right attributes and apply normalization if it's active
    for (i = x->x_attributeLo, j = 0; i <= x->x_attributeHi; i++, j++)
    {
        t_attributeIdx thisAttribute;
        thisAttribute = x->x_attributeData[i].order;

        if (thisAttribute > instance1.length - 1 || thisAttribute > instance2.length - 1)
        {
            pd_error (x, "%s: attribute %i does not exist for both feature vectors. cannot compute distance.", x->x_objSymbol->s_name, thisAttribute);

            // free local memory before exiting
            t_freebytes (vec1Buffer, vecLen * sizeof (t_float));
            t_freebytes (vec2Buffer, vecLen * sizeof (t_float));
            t_freebytes (vecWeights, vecLen * sizeof (t_float));

            return (FLT_MAX);
        }

        vecWeights[j] = x->x_attributeData[thisAttribute].weight;

        if (x->x_normalize)
        {
            // don't need to test if instance1 has attribute values below/above current min/max values in normData, because this function isn't used on input vectors from the wild. It's only used to compare database instances to other database instances, so the normalization terms are valid.
            min = x->x_attributeData[thisAttribute].normData.minVal;
            max = x->x_attributeData[thisAttribute].normData.maxVal;
            normScalar = x->x_attributeData[thisAttribute].normData.normScalar;

            vec1Buffer[j] = (instance1.data[thisAttribute] - min) * normScalar;
            vec2Buffer[j] = (instance2.data[thisAttribute] - min) * normScalar;
        }
        else
        {
            vec1Buffer[j] = instance1.data[thisAttribute];
            vec2Buffer[j] = instance2.data[thisAttribute];
        }
    }

    switch (x->x_distMetric)
    {
        case euclidean:
            dist = tIDLib_euclidDist (vecLen, vec1Buffer, vec2Buffer, vecWeights, false);
            break;
        case taxi:
            dist = tIDLib_taxiDist (vecLen, vec1Buffer, vec2Buffer, vecWeights);
            break;
        case correlation:
            dist = tIDLib_corr (vecLen, vec1Buffer, vec2Buffer);
            // bash to the 0-2 range, then flip sign so that lower is better. this keeps things consistent with other distance metrics.
            dist += 1;
            dist *= -1;
            break;
        default:
            break;
    };

    // free local memory
    t_freebytes (vec1Buffer, vecLen * sizeof (t_float));
    t_freebytes (vec2Buffer, vecLen * sizeof (t_float));
    t_freebytes (vecWeights, vecLen * sizeof (t_float));

    return (dist);
}


// this gives the distance between an input feature vector and a defined class reference
static t_float timbreID_getClassRefDist (t_timbreID* x, t_attributeIdx inputLen, t_float* input, t_instance classRef)
{
    t_attributeIdx i, classLen;
    t_float dist;
    t_float* weights;

    classLen = classRef.length;
    dist = 0.0;

    weights = (t_float *)t_getbytes (classLen * sizeof (t_float));

    for (i = 0; i < classLen; i++)
    {
        if (i >= inputLen)
        {
            pd_error (x, "%s: attribute %i does not exist for both feature vectors. cannot compute distance.", x->x_objSymbol->s_name, i);

            // free local memory before exiting
            t_freebytes (weights, classLen * sizeof (t_float));

            // TODO: shouldn't have an else branch after a return statement
            return (FLT_MAX);
        }
        else
            weights[i] = 1.0;
    }

    switch (x->x_distMetric)
    {
        case euclidean:
            dist = tIDLib_euclidDist (classLen, input, classRef.data, weights, false);
            break;
        case taxi:
            dist = tIDLib_taxiDist (classLen, input, classRef.data, weights);
            break;
        case correlation:
            dist = tIDLib_corr (classLen, input, classRef.data);
            // bash to the 0-2 range, then flip sign so that lower is better. this keeps things consistent with other distance metrics.
            dist += 1;
            dist *= -1;
            break;
        default:
            break;
    };

    // free local memory
    t_freebytes (weights, classLen * sizeof (t_float));

    return (dist);
}


static void timbreID_attributeDataResize (t_timbreID* x, t_attributeIdx oldSize, t_attributeIdx newSize, t_bool postFlag)
{
    t_attributeIdx i;

    x->x_attributeData = (t_attributeData *)t_resizebytes (x->x_attributeData, oldSize * sizeof (t_attributeData), newSize * sizeof (t_attributeData));
    x->x_maxFeatureLength = newSize;

    if (x->x_maxFeatureLength < x->x_minFeatureLength)
        x->x_minFeatureLength = x->x_maxFeatureLength;

    // initialize attributeData
    for (i = 0; i < x->x_maxFeatureLength; i++)
    {
        x->x_attributeData[i].inputData = 0.0;
        x->x_attributeData[i].order = i;
        x->x_attributeData[i].weight = 1.0;
        x->x_attributeData[i].name = gensym ("NULL");
    }

    x->x_normalize = false;
    x->x_attributeLo = 0;
    x->x_attributeHi = x->x_maxFeatureLength - 1;

    if (postFlag)
    {
        post ("%s: max feature length: %i.", x->x_objSymbol->s_name, x->x_maxFeatureLength);
        post ("%s: min feature length: %i.", x->x_objSymbol->s_name, x->x_minFeatureLength);
        post ("%s: feature attribute normalization OFF.", x->x_objSymbol->s_name);
        post ("%s: attribute weights initialized", x->x_objSymbol->s_name);
        post ("%s: attribute order initialized", x->x_objSymbol->s_name);
        post ("%s: attribute range: %i through %i.", x->x_objSymbol->s_name, x->x_attributeLo, x->x_attributeHi);
    }
};
/* -------- END utility functions -------- */




/* -------- classification methods -------- */
static void timbreID_train (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    t_instanceIdx i, instanceIdx;

    if (x->x_normalize)
        pd_error (x, "%s: cannot add more training instances when database is normalized. deactivate normalization first.", x->x_objSymbol->s_name);
    else if (x->x_numClusters != x->x_numInstances)
        pd_error (x, "%s: cannot add more training instances when database is clustered. uncluster first.", x->x_objSymbol->s_name);
    else
    {
    instanceIdx = x->x_numInstances;

    x->x_instances = (t_instance *)t_resizebytes (x->x_instances, x->x_numInstances * sizeof (t_instance), (x->x_numInstances + 1) * sizeof (t_instance));
    x->x_clusters = (t_cluster *)t_resizebytes (x->x_clusters, x->x_numInstances * sizeof (t_cluster), (x->x_numInstances + 1) * sizeof (t_cluster));

    x->x_instances[instanceIdx].clusterMembership = instanceIdx;
    x->x_instances[instanceIdx].length = argc;

    x->x_instances[instanceIdx].data = (t_float *)t_getbytes (x->x_instances[instanceIdx].length * sizeof (t_float));

    x->x_clusters[instanceIdx].numMembers = 2; // 2 because we're unclustered to start, and each instance has a cluster with itself as a member, plus UINT_MAX as the 2nd element to terminate the list

    x->x_clusters[instanceIdx].members = (t_instanceIdx *)t_getbytes (x->x_clusters[instanceIdx].numMembers * sizeof (t_instanceIdx));

    // init new clusterMembers
    x->x_clusters[instanceIdx].members[0] = instanceIdx; // first member of the cluster is the instance index
    x->x_clusters[instanceIdx].members[1] = UINT_MAX;

    x->x_clusters[instanceIdx].votes = 0;

    x->x_numInstances++;
    x->x_numClusters++;
    x->x_neighborhood++;

    if (x->x_instances[instanceIdx].length > x->x_maxFeatureLength)
        timbreID_attributeDataResize (x, x->x_maxFeatureLength, x->x_instances[instanceIdx].length, true);

    if (x->x_instances[instanceIdx].length < x->x_minFeatureLength)
        x->x_minFeatureLength = x->x_instances[instanceIdx].length;

    for (i = 0; i < x->x_instances[instanceIdx].length; i++)
        x->x_instances[instanceIdx].data[i] = atom_getfloat (argv + i);

    outlet_float (x->x_id, instanceIdx); // output "received" feedback here rather than a post (in case of hi-speed training)
    }
}


static void timbreID_id (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    t_float dist, bestDist, secondBestDist, confidence;
    t_instanceIdx i, winningID, topVoteInstances;
    t_uInt topVote;
    t_attributeIdx listLength;

    if (x->x_numInstances)
    {
        listLength = argc;
        confidence = 0.0;

        if (listLength > x->x_maxFeatureLength)
        {
            pd_error (x, "%s: input feature list longer than current max feature length of database. input ignored.", x->x_objSymbol->s_name);
            return;
        };

        // init votes to 0
        for (i = 0; i < x->x_numClusters; i++)
            x->x_clusters[i].votes = 0;

        // init cluster info to instance idx
        for (i = 0; i < x->x_numInstances; i++)
        {
            x->x_instances[i].knnInfo.idx = i;
            x->x_instances[i].knnInfo.cluster = i;
            x->x_instances[i].knnInfo.dist = FLT_MAX;
            x->x_instances[i].knnInfo.safeDist = FLT_MAX;
        }

        for (i = 0; i < listLength; i++)
            x->x_attributeData[i].inputData = atom_getfloat (argv + i);

        winningID = UINT_MAX;
        bestDist = FLT_MAX;
        secondBestDist = FLT_MAX;

        for (i = 0; i < x->x_numInstances; i++)
        {
            dist = timbreID_getInputDist (x, i);

            // abort _id() altogether if distance measurement fails. return value of FLT_MAX indicates failure
            if (dist == FLT_MAX)
                return;

            x->x_instances[i].knnInfo.dist = x->x_instances[i].knnInfo.safeDist = dist; // store the distance
        };

        // a reduced sort, so that the first k elements in x_instances will be the ones with lowest distances to the input feature in the list, and in order to boot. NOTE: tIDLib_sortKnnInfo() does not operate on the primary instance data, just each instance's knnInfo member.
        tIDLib_sortKnnInfo (x->x_k, x->x_numInstances, UINT_MAX, x->x_instances); // pass a prevMatch value of UINT_MAX as a flag that it's unused here

        // store instance's cluster id
        for (i = 0; i < x->x_k; i++)
        {
            t_instanceIdx thisInstance;
            thisInstance = x->x_instances[i].knnInfo.idx;

            x->x_instances[i].knnInfo.cluster = x->x_instances[thisInstance].clusterMembership;
        }

        // vote
        for (i = 0; i < x->x_k; i++)
        {
            t_instanceIdx thisCluster;
            thisCluster = x->x_instances[i].knnInfo.cluster;

            x->x_clusters[thisCluster].votes++;
        }

        // TODO: potential issue is that i'm using t_uInt for t_cluster.votes and topVote, so can no longer have -1 as the initialized value. should be that initializing to 0 is actually ok, even in the case of a tie for 0. originally this was topVote = -1;
        topVote = 0;

        for (i = 0; i < x->x_numClusters; i++)
        {
            if (x->x_clusters[i].votes > topVote)
            {
                topVote = x->x_clusters[i].votes;
                winningID = i; // store cluster id of winner
            };
        }

        // see if there is a tie for number of votes
        topVoteInstances = 0;

        for (i = 0; i < x->x_numClusters; i++)
        {
            if (x->x_clusters[i].votes == topVote)
                topVoteInstances++;
        }

        // in case of a tie, pick the instance with the shortest distance. The knnInfo for instances will be sorted by distance, so x_instance[0].knnInfo.cluster is the cluster ID of the instance with the smallest distance
        if (topVoteInstances > 1)
            winningID = x->x_instances[0].knnInfo.cluster;

        for (i = 0; i < x->x_k; i++)
        {
            if (x->x_instances[i].knnInfo.cluster == winningID)
            {
                bestDist = x->x_instances[i].knnInfo.safeDist;
                break;
            };
        }

        // this assigns the distance belonging to the first item in knnInfo that isn't a member of the winning cluster to the variable "secondBest". i.e., the distance between the test vector and the next nearest instance that isn't a member of the winning cluster.
        for (i = 0; i < x->x_k; i++)
        {
            if (x->x_instances[i].knnInfo.cluster != winningID)
            {
                secondBestDist = x->x_instances[i].knnInfo.safeDist;
                break;
            };
        }

    /*
    // unless I'm missing something, this is faulty logic with the current setup. If all K items belong to the same cluster, then secondBestDist will still be FLT_MAX, which will make confidence below near to 1.0. That's good. Assigning it the actual second best distance from a sibling cluster member will produce a number near 0 for confidence, which is wrong.

        // if no second best assignment is made (because all K items belong to same cluster), make 2nd best the 2nd in list
        if (secondBestdist == FLT_MAX)
            secondBestDist = x->x_instances[1].knnInfo.safeDist;

    // if above note on faulty logic is correct, then this is not needed. confidence can just be calculated without a condition
        if (secondBestDist <= 0)
            confidence = -FLT_MAX;
        else
            confidence = 1.0-(bestDist/secondBestDist);
    */

        confidence = 1.0 - (bestDist / secondBestDist);

        outlet_float (x->x_confidence, confidence);
        outlet_float (x->x_nearestDist, bestDist);
        outlet_float (x->x_id, winningID);

        if (x->x_outputKnnMatches)
        {
            for (i = 0; i < x->x_k; i++)
            {
                // suppress reporting of the vote-winning id, because it was already reported above
                if (x->x_instances[i].knnInfo.cluster != winningID)
                {
                    outlet_float (x->x_nearestDist, x->x_instances[i].knnInfo.safeDist);
                    outlet_float (x->x_id, x->x_instances[i].knnInfo.cluster);
                }
            }
        }
    }
    else
        pd_error (x, "%s: no training instances have been loaded. cannot perform ID.", x->x_objSymbol->s_name);
}


static void timbreID_worstMatch (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    t_float dist, worstDist;
    t_attributeIdx listLength;
    t_instanceIdx i, losingID;

    if (x->x_numInstances)
    {
        listLength = argc;

        if (listLength > x->x_maxFeatureLength)
        {
            pd_error (x, "%s: input feature list too long. input ignored.", x->x_objSymbol->s_name);
            return;
        };

        for (i = 0; i < listLength; i++)
            x->x_attributeData[i].inputData = atom_getfloat (argv + i);

        losingID = UINT_MAX;
        worstDist = -FLT_MAX;

        for (i = 0; i < x->x_numInstances; i++)
        {
            dist = timbreID_getInputDist (x, i);

            // abort _worstMatch() altogether if distance measurement fails
            if (dist == FLT_MAX)
                return;

            if (dist > worstDist)
            {
                worstDist = dist;
                losingID = i;
            }
        };

        losingID = x->x_instances[losingID].clusterMembership;

        outlet_float (x->x_confidence, 0);
        outlet_float (x->x_nearestDist, worstDist);
        outlet_float (x->x_id, losingID);
    }
    else
        pd_error (x, "%s: no training instances have been loaded. cannot perform worst match.", x->x_objSymbol->s_name);
}
/* -------- END classification methods -------- */




/* -------- concatenative synthesis methods -------- */
static void timbreID_concatId (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    if (x->x_numInstances)
    {
        t_float dist, bestDist;
        t_instanceIdx i, j, halfNeighborhood, winningID;
        t_sLongInt searchStart; // need this to be signed for wraparound
        t_attributeIdx listLength;

        // init knn info
        for (i = 0; i < x->x_numInstances; i++)
        {
            x->x_instances[i].knnInfo.idx = i;
            x->x_instances[i].knnInfo.dist = x->x_instances[i].knnInfo.safeDist = FLT_MAX;
        }

        listLength = argc;

        halfNeighborhood = x->x_neighborhood * 0.5;

        if (listLength > x->x_maxFeatureLength)
        {
            pd_error (x, "%s: input feature list too long. input ignored.", x->x_objSymbol->s_name);
            return;
        };

        for (i = 0; i < listLength; i++)
            x->x_attributeData[i].inputData = atom_getfloat (argv + i);

        winningID = UINT_MAX;
        bestDist = FLT_MAX;

        // searchStart needs to be signed in order to do wraparound. since t_instanceIdx is a unsigned int, we need a signed long int. unsigned int overflow is undefined, so can't count on that as an alternative for wraparound into negative
        if (x->x_concatWrap)
        {
            // for just the searchStart, check to see if x->x_neighborhood was EVEN.  if so, we should make searchStart = searchCenter - halfNeighborhood + 1
            if (x->x_neighborhood % 2 == 0)
                searchStart = x->x_searchCenter - halfNeighborhood + 1;
            else
                searchStart = x->x_searchCenter - halfNeighborhood;

            if (searchStart < 0)
                searchStart = x->x_numInstances + searchStart;  // wraps in reverse to end of table.  + searchStart because searchStart is negative here.
            else
                searchStart = searchStart % x->x_numInstances;
        }
        else
        {
            // for just the searchStart, check to see if x->x_neighborhood was EVEN.  if so, we should make searchStart = searchCenter - halfNeighborhood + 1
            if (x->x_neighborhood % 2 == 0)
                searchStart = x->x_searchCenter - halfNeighborhood + 1;
            else
                searchStart = x->x_searchCenter - halfNeighborhood;

            // TODO: is this (t_instanceIdx) typecast a good idea here? Just to get rid of warning?
            if (searchStart < 0)
                searchStart = 0; // no wrapping
            else if ((t_instanceIdx)searchStart >= x->x_numInstances)
                searchStart = x->x_numInstances - 1;
        }

        for (j = 0, i = searchStart; j < x->x_neighborhood; j++)
        {
            dist = timbreID_getInputDist (x, i);

            // abort _concat_id() altogether if distance measurement fails
            if (dist == FLT_MAX)
                return;

            x->x_instances[i].knnInfo.dist = x->x_instances[i].knnInfo.safeDist = dist; // store the distance

             i++;

             if (x->x_concatWrap)
                 i = i % x->x_numInstances;
             else
             {
                 if (i >= x->x_numInstances)
                     break;
             }
        };

        // a reduced sort, so that the first maxMatches elements in x_instances will have the lowest distances in their knnInfo.safeDist, and in order to boot.
        // pass x->x_prevMatch to make sure we don't output the same match two times in a row (to prevent one grain being played back several
        // times in sequence.
        // this is wasteful in restricted searches because we don't need to look through all x->x_numInstances
        tIDLib_sortKnnInfo (x->x_maxMatches, x->x_numInstances, x->x_prevMatch, x->x_instances);

        if (x->x_prevMatch == UINT_MAX)
        {
            winningID = x->x_instances[0].knnInfo.idx;
            bestDist = x->x_instances[0].knnInfo.safeDist;
        }
        else
        {
            for (i = 0, bestDist = FLT_MAX; i < x->x_maxMatches; i++)
            {
                t_instanceIdx thisInstance;

                thisInstance = x->x_instances[i].knnInfo.idx;

                dist = timbreID_getDist (x, x->x_instances[x->x_prevMatch], x->x_instances[thisInstance]);

                // abort _concat_id() altogether if distance measurement fails
                if (dist == FLT_MAX)
                    return;

                if (dist < bestDist)
                {
                    bestDist = dist;
                    winningID = thisInstance;
                };
            };
        };

        if (x->x_reorientFlag)
            x->x_searchCenter = winningID;

        if (rand() < RAND_MAX * x->x_jumpProb)
        {
            winningID += rand();
            winningID = winningID % x->x_numInstances;
            x->x_searchCenter = winningID;
        };

        if (x->x_stutterProtect)
            x->x_prevMatch = winningID;
        else
            x->x_prevMatch = UINT_MAX;

        outlet_float (x->x_nearestDist, bestDist);
        outlet_float (x->x_id, winningID);
    }
    else
        pd_error (x, "%s: no training instances have been loaded. cannot perform ID.", x->x_objSymbol->s_name);
}


static void timbreID_concatSearchWrap (t_timbreID* x, t_floatarg flag)
{
    flag = (flag > 1) ? 1 : flag;
    flag = (flag < 0) ? 0 : flag;
    x->x_concatWrap = flag;
}


static void timbreID_concatNeighborhood (t_timbreID* x, t_floatarg n)
{
    n = (n > x->x_numInstances) ? x->x_numInstances : n;
    n = (n < 1) ? 1 : n;
    x->x_neighborhood = n;
}


static void timbreID_concatJumpProb (t_timbreID* x, t_floatarg jp)
{
    jp = (jp < 0) ? 0 : jp;
    jp = (jp > 1) ? 1 : jp;
    x->x_jumpProb = jp;
}


static void timbreID_concatReorient (t_timbreID* x, t_floatarg r)
{
    r = (r < 0) ? 0 : r;
    r = (r > 1) ? 1 : r;
    x->x_reorientFlag = r;
}


static void timbreID_concatSearchCenter (t_timbreID* x, t_floatarg sc)
{
    if (sc < 0)
        x->x_searchCenter = 0;
    else if (sc >= x->x_numInstances)
        x->x_searchCenter = x->x_numInstances - 1;
    else
        x->x_searchCenter = sc;
}


static void timbreID_concatMaxMatches (t_timbreID* x, t_floatarg mm)
{
    // do this safety check on the incoming floatarg, since doing it later on x_maxMatches won't work. x_maxMatches is an unsigned long, so can't handle negative values normally
    mm = (mm < 1) ? 1 : mm;

    x->x_maxMatches = (mm > x->x_numInstances) ? x->x_numInstances : mm;
// x->x_maxMatches = (x->x_maxMatches<1)?1:x->x_maxMatches;
}


static void timbreID_concatStutterProtect (t_timbreID* x, t_floatarg sp)
{
    sp = (sp < 0) ? 0 : sp;
    sp = (sp > 1) ? 1 : sp;
    x->x_stutterProtect = sp;
}
/* -------- END concatenative synthesis methods -------- */




static void timbreID_knn (t_timbreID* x, t_floatarg k)
{
    if (k < 1.0)
        pd_error (x, "%s: k must be greater than zero.", x->x_objSymbol->s_name);
    else if (k > x->x_numInstances)
        pd_error (x, "%s: k must be less than the total number of instances.", x->x_objSymbol->s_name);
    else
    {
        x->x_k = k;
        post ("%s: searching %i neighbors for KNN.", x->x_objSymbol->s_name, x->x_k);
    }
}


static void timbreID_outputKnnMatches (t_timbreID* x, t_floatarg out)
{
    out = (out < 0) ? 0 : out;
    out = (out > 1) ? 1 : out;

    x->x_outputKnnMatches = out;

    if (x->x_outputKnnMatches)
        post ("%s: reporting all KNN matches.", x->x_objSymbol->s_name);
    else
        post ("%s: reporting best match only.", x->x_objSymbol->s_name);
}


static void timbreID_normalize (t_timbreID* x, t_floatarg n)
{
    t_float* attributeColumn;

    // create local memory
    attributeColumn = (t_float *)t_getbytes (x->x_numInstances * sizeof (t_float));

    if (n <= 0)
    {
        t_attributeIdx i;

        // initialize normData
        for (i = 0; i < x->x_maxFeatureLength; i++)
        {
            x->x_attributeData[i].normData.minVal = 0.0;
            x->x_attributeData[i].normData.maxVal = 0.0;
            x->x_attributeData[i].normData.normScalar = 1.0;
        }

        x->x_normalize = false;
        post ("%s: feature attribute normalization OFF.", x->x_objSymbol->s_name);
    }
    else
    {
        if (x->x_numInstances)
        {
            t_instanceIdx i;
            t_attributeIdx j;

            // j for columns (attributes), i for rows (instances)
            for (j = 0; j < x->x_maxFeatureLength; j++)
            {
                for (i = 0; i < x->x_numInstances; i++)
                {
                    if (j > x->x_instances[i].length - 1)
                    {
                        pd_error (x, "%s: attribute %i out of range for database instance %i. aborting normalization", x->x_objSymbol->s_name, j, i);
                        // initialize normData
                        for (i = 0; i < x->x_maxFeatureLength; i++)
                        {
                            x->x_attributeData[i].normData.minVal = 0.0;
                            x->x_attributeData[i].normData.maxVal = 0.0;
                            x->x_attributeData[i].normData.normScalar = 1.0;
                        }

                        x->x_normalize = false;

                        post ("%s: feature attribute normalization OFF.", x->x_objSymbol->s_name);
                        // free local memory before exit
                        t_freebytes (attributeColumn, x->x_numInstances * sizeof (t_float));
                        // TODO: shouldn't have an else branch after a return statement
                        return;
                    }
                    else
                        attributeColumn[i] = x->x_instances[i].data[j];
                }

                tIDLib_bubbleSort (x->x_numInstances, attributeColumn);

                x->x_attributeData[j].normData.minVal = attributeColumn[0];
                x->x_attributeData[j].normData.maxVal = attributeColumn[x->x_numInstances - 1];

                if (x->x_attributeData[j].normData.maxVal <= x->x_attributeData[j].normData.minVal)
                {
                    // this will fix things in the case of 1 instance, where min==max
                    x->x_attributeData[j].normData.minVal = 0.0;
                    x->x_attributeData[j].normData.normScalar = 1.0;
                }
                else
                    x->x_attributeData[j].normData.normScalar = 1.0 / (x->x_attributeData[j].normData.maxVal - x->x_attributeData[j].normData.minVal);
            };

            x->x_normalize = true;
            post ("%s: feature attribute normalization ON.", x->x_objSymbol->s_name);
        }
        else
            pd_error (x, "%s: no training instances have been loaded. cannot calculate normalization terms.", x->x_objSymbol->s_name);
    }

    // free local memory
    t_freebytes (attributeColumn, x->x_numInstances * sizeof (t_float));

}


static void timbreID_print (t_timbreID* x)
{
    post ("timbreID version %s", TID_VERSION);
    post ("%s: no. of instances: %i", x->x_objSymbol->s_name, x->x_numInstances);
    post ("%s: max feature length: %i", x->x_objSymbol->s_name, x->x_maxFeatureLength);
    post ("%s: min feature length: %i", x->x_objSymbol->s_name, x->x_minFeatureLength);
    post ("%s: attribute range: %i through %i", x->x_objSymbol->s_name, x->x_attributeLo, x->x_attributeHi);
    post ("%s: normalization: %i", x->x_objSymbol->s_name, x->x_normalize);
    post ("%s: distance metric: %i", x->x_objSymbol->s_name, x->x_distMetric);
    post ("%s: no. of clusters: %i", x->x_objSymbol->s_name, x->x_numClusters);
    post ("%s: no. of defined reference classes: %i", x->x_objSymbol->s_name, x->x_numClassRefs);
    post ("%s: KNN: %i", x->x_objSymbol->s_name, x->x_k);
    post ("%s: output KNN matches: %i", x->x_objSymbol->s_name, x->x_outputKnnMatches);
    post ("%s: relative ordering: %i", x->x_objSymbol->s_name, x->x_relativeOrdering);
    post ("%s: concatenative wraparound: %i", x->x_objSymbol->s_name, x->x_concatWrap);
    post ("%s: search center: %i", x->x_objSymbol->s_name, x->x_searchCenter);
    post ("%s: neighborhood: %i", x->x_objSymbol->s_name, x->x_neighborhood);
    post ("%s: reorient: %i", x->x_objSymbol->s_name, x->x_reorientFlag);
    post ("%s: max matches: %i", x->x_objSymbol->s_name, x->x_maxMatches);
    post ("%s: jump probability: %f", x->x_objSymbol->s_name, x->x_jumpProb);
    post ("%s: stutter protect: %i", x->x_objSymbol->s_name, x->x_stutterProtect);
}


static void timbreID_numInstances (t_timbreID* x)
{
    t_atom listOut;
    t_symbol* selector;

    selector = gensym ("num_instances");

    SETFLOAT (&listOut, x->x_numInstances);

    outlet_anything (x->x_listOut, selector, 1, &listOut);
}


static void timbreID_manualCluster (t_timbreID* x, t_floatarg numClusters, t_floatarg clusterIdx, t_floatarg low, t_floatarg hi)
{
    t_instanceIdx i, j, clusterIdxInt, lowIdx, hiIdx, numMembers;

    clusterIdxInt = (clusterIdx < 0) ? 0 : clusterIdx;
    lowIdx = (low < 0) ? 0 : low;
    hiIdx = (hi < 0) ? 0 : hi;

    if (lowIdx > hiIdx)
    {
        t_instanceIdx tmp;

        tmp = lowIdx;
        lowIdx = hiIdx;
        hiIdx = tmp;
    }

    numMembers = hiIdx - lowIdx + 1;

    if (lowIdx > x->x_numInstances - 1 || hiIdx > x->x_numInstances - 1)
    {
        pd_error (x, "%s: instances out of range. clustering failed.", x->x_objSymbol->s_name);
        return;
    }

    if (clusterIdx >= numClusters)
    {
        pd_error (x, "%s: cluster index %i out of range for given number of clusters %i. clustering failed.", x->x_objSymbol->s_name, clusterIdxInt, (t_instanceIdx)numClusters);
        return;
    }

    if (x->x_numInstances < numClusters)
        pd_error (x, "%s: not enough instances to cluster. clustering failed.", x->x_objSymbol->s_name);
    else
    {
        // update the number of clusters, trusting that the user knows what they're doing!
        x->x_numClusters = numClusters;

        // free this cluster's members memory
        t_freebytes (x->x_clusters[clusterIdxInt].members, x->x_clusters[clusterIdxInt].numMembers * sizeof (t_instanceIdx));

        // update the .numMembers value
        x->x_clusters[clusterIdxInt].numMembers = numMembers + 1; // +1 for the terminating UINT_MAX

        // get new memory for this cluster's members
        x->x_clusters[clusterIdxInt].members = (t_instanceIdx *)t_getbytes (x->x_clusters[clusterIdxInt].numMembers * sizeof (t_instanceIdx));

        x->x_clusters[clusterIdxInt].votes = 0;

        for (i = lowIdx; i <= hiIdx; i++)
            x->x_instances[i].clusterMembership = clusterIdxInt;

        for (i = lowIdx, j = 0; i <= hiIdx; i++, j++)
            x->x_clusters[clusterIdxInt].members[j] = i;

        // terminate with UINT_MAX
        x->x_clusters[clusterIdxInt].members[j] = UINT_MAX;

        // resize the excess clusterMembers memory back to the default size of 2 and store the default instance index as the cluster member followed by a terminating UINT_MAX
        for (i = x->x_numClusters; i < x->x_numInstances; i++)
        {
            x->x_clusters[i].members = (t_instanceIdx *)t_resizebytes (x->x_clusters[i].members, x->x_clusters[i].numMembers * sizeof (t_instanceIdx), 2 * sizeof (t_instanceIdx));
            x->x_clusters[i].numMembers = 2;
            x->x_clusters[i].members[0] = i;
            x->x_clusters[i].members[1] = UINT_MAX;
        }

        post ("%s: cluster %i contains instances %i through %i.", x->x_objSymbol->s_name, clusterIdxInt, lowIdx, hiIdx);
    };
}


static void timbreID_computeCluster (t_timbreID* x, t_floatarg numClusters)
{
    if (x->x_numInstances < numClusters)
        pd_error (x, "%s: not enough instances to cluster.", x->x_objSymbol->s_name);
    else if (x->x_numClusters != x->x_numInstances)
        pd_error (x, "%s: instances already clustered. uncluster first.", x->x_objSymbol->s_name);
    else if (numClusters == 0)
        pd_error (x, "%s: cannot create 0 clusters.", x->x_objSymbol->s_name);
    else
    {
    t_instanceIdx i, j, k, numInstances, numInstancesM1, numPairs, numClusterMembers1, numClusterMembers2, numClusterMembersSum, clusterCount;
    t_instanceIdx* minDistIdx;

    t_instance* clusterData;
    t_float* pairDists;
    t_float minDist, numClusterMembers1_recip;
    t_atom* listOut;
    t_symbol* selector;

    x->x_numClusters = numClusters;
    numInstances = x->x_numInstances;
    numInstancesM1 = numInstances - 1;
    numPairs = (numInstances * numInstancesM1) * 0.5;
    clusterCount = numInstances;
    numClusterMembers1 = 0;
    numClusterMembers2 = 0;
    numClusterMembers1_recip = 1;
    i = j = k = 0;

    // create local memory
    minDistIdx = (t_instanceIdx *)t_getbytes (2 * sizeof (t_instanceIdx));
    clusterData = (t_instance *)t_getbytes (numInstances * sizeof (t_instance));
    pairDists = (t_float *)t_getbytes (numPairs * sizeof (t_float));
    listOut = (t_atom *)t_getbytes (numInstances * sizeof (t_atom));

    for (i = 0; i < numInstances; i++)
    {
        x->x_clusters[i].members[0] = i; // first member of the cluster is the instance index
        x->x_clusters[i].members[1] = UINT_MAX;
    }

    // copy x->x_instances into a safe local copy: clusterData
    for (i = 0; i < numInstances; i++)
    {
        clusterData[i].length = x->x_instances[i].length;
        clusterData[i].data = (t_float *)t_getbytes (x->x_instances[i].length * sizeof (t_float));

        for (j = 0; j < clusterData[i].length; j++)
            clusterData[i].data[j] = x->x_instances[i].data[j];
    }

    while (clusterCount > x->x_numClusters)
    {
        minDist = FLT_MAX;

        // init minDistIdx to UINT_MAX as an indicator that it's initialized
        for (i = 0; i < 2; i++)
            minDistIdx[i] = UINT_MAX;

        // init pair distances
        for (i = 0; i < numPairs; i++)
            pairDists[i] = FLT_MAX;

        // get distances between all possible pairs in clusterData
        for (i = 0, k = 0; i < numInstancesM1; i++)
        {
            if (clusterData[i].data[0] != FLT_MAX) // if this is true, the data hasn't been clustered yet.
            {
                for (j = 1; j < numInstances; j++)
                {
                    if ((i + j) < numInstances)
                    {
                        if (clusterData[i+j].data[0] != FLT_MAX)
                        {
                            pairDists[k] = timbreID_getDist (x, clusterData[i], clusterData[i + j]);
                            numClusterMembers1 = x->x_clusters[i].numMembers - 1; // -1 because the list is terminated with UINT_MAX
                            numClusterMembers2 = x->x_clusters[i+j].numMembers - 1;

                            // definition of Ward's linkage from MATLAB linkage doc
                            // pairDists[k] is already euclidean distance
                            numClusterMembersSum = numClusterMembers1 + numClusterMembers2;

                            if (numClusterMembersSum > 0)
                                pairDists[k] = numClusterMembers1 * numClusterMembers2 * (pairDists[k] / (numClusterMembers1+numClusterMembers2));
                            else
                                pairDists[k] = FLT_MAX;

                            if (pairDists[k] < minDist)
                            {
                                minDist = pairDists[k];
                                minDistIdx[0] = i;
                                minDistIdx[1] = i + j;
                            };

                            k++; // increment pairDists index if something was actually written to it.
                        };
                    }
                    else
                        break;
                }
            }
        };

        // we've found the smallest distance between clusterData elements and stored it
        // in minDist. we've store the clusterData indices of the two elements in
        // minDistIdx[0] and minDistIdx[1].

        // set i to the index for storing the new member (s) of the cluster.
        i = x->x_clusters[minDistIdx[0]].numMembers - 1;

        // actually store the new member (s).
        j = 0;

        while (x->x_clusters[minDistIdx[1]].members[j] != UINT_MAX)
        {
            // make some more memory for the new member (s)
            x->x_clusters[minDistIdx[0]].members = (t_instanceIdx *)t_resizebytes (x->x_clusters[minDistIdx[0]].members, x->x_clusters[minDistIdx[0]].numMembers * sizeof (t_instanceIdx), (x->x_clusters[minDistIdx[0]].numMembers + 1) * sizeof (t_instanceIdx));
            (x->x_clusters[minDistIdx[0]].numMembers)++; // remember to update this member list's length

            x->x_clusters[minDistIdx[0]].members[i++] = x->x_clusters[minDistIdx[1]].members[j++];
        }

        i = x->x_clusters[minDistIdx[0]].numMembers - 1;
        x->x_clusters[minDistIdx[0]].members[i] = UINT_MAX; // terminate with UINT_MAX

        numClusterMembers1 = x->x_clusters[minDistIdx[0]].numMembers - 1;

        if (numClusterMembers1 > 0)
            numClusterMembers1_recip = 1.0 / (t_float)numClusterMembers1;
        else
            numClusterMembers1_recip = 1.0;

        // resize the usurped cluster's cluster list memory, and update its size to 1
        x->x_clusters[minDistIdx[1]].members = (t_instanceIdx *)t_resizebytes (x->x_clusters[minDistIdx[1]].members, x->x_clusters[minDistIdx[1]].numMembers * sizeof (t_instanceIdx), sizeof (t_instanceIdx));
        x->x_clusters[minDistIdx[1]].members[0] = UINT_MAX; // terminate with UINT_MAX
        x->x_clusters[minDistIdx[1]].numMembers = 1;

        // grab the first original instance for this cluster index
        for (i = 0; i < clusterData[minDistIdx[0]].length; i++)
            clusterData[minDistIdx[0]].data[i] = x->x_instances[minDistIdx[0]].data[i];

        // sum the original instances of the cluster members to compute centroid below
        for (i = 1; i < numClusterMembers1; i++)
            for (j = 0; j < clusterData[minDistIdx[0]].length; j++)
                clusterData[minDistIdx[0]].data[j] += x->x_instances[  x->x_clusters[minDistIdx[0]].members[i]  ].data[j];

        // compute centroid
        for (i = 0; i < clusterData[minDistIdx[0]].length; i++)
            clusterData[minDistIdx[0]].data[i] *= numClusterMembers1_recip;

        // write FLT_MAX to the first element in the nearest neighbor's instance to indicate it's now vacant.
        // this is all that's needed since all previous members were averaged and stored here.
        clusterData[minDistIdx[1]].data[0] = FLT_MAX;

        clusterCount--;
    };

    // since the indices of the clusters have gaps from the process,
    // shift the clusterMembers arrays that actually have content ( != UINT_MAX)
    // to the head of clusterMembers.  this will produce indices from 0 through numClusters - 1.
    for (i = 0, k = 0; i < numInstances; i++)
    {
        if (x->x_clusters[i].members[0] != UINT_MAX)
        {
            // resize this member list
             x->x_clusters[k].members = (t_instanceIdx *)t_resizebytes (x->x_clusters[k].members, x->x_clusters[k].numMembers * sizeof (t_instanceIdx), x->x_clusters[i].numMembers * sizeof (t_instanceIdx));

            for (j = 0; j < x->x_clusters[i].numMembers; j++)
                x->x_clusters[k].members[j] = x->x_clusters[i].members[j];

            // shift the list length info back
             x->x_clusters[k].numMembers = x->x_clusters[i].numMembers;

            k++;
        };
    }

    // resize the excess clusterMembers memory back to the default size of 2 and store the default instance index as the cluster member followed by a terminating UINT_MAX
    for (i = x->x_numClusters; i < numInstances; i++)
    {
        x->x_clusters[i].members = (t_instanceIdx *)t_resizebytes (x->x_clusters[i].members, x->x_clusters[i].numMembers * sizeof (t_instanceIdx), 2 * sizeof (t_instanceIdx));
        x->x_clusters[i].numMembers = 2;
        x->x_clusters[i].members[0] = i;
        x->x_clusters[i].members[1] = UINT_MAX;
    }

    // UPDATE: as of version 0.7, not doing this resize because we're always keeping x_clusters to be size x_numInstances
    // resize clusterMembers so it is only x->x_numClusters big
// x->x_clusters = (t_cluster *)t_resizebytes (x->x_clusters, numInstances * sizeof (t_cluster), x->x_numClusters * sizeof (t_cluster));

    for (i = 0, k = 0; i < x->x_numClusters; i++)
        for (j = 0; j < x->x_clusters[i].numMembers - 1; j++, k++)
        {
            x->x_instances[x->x_clusters[i].members[j]].clusterMembership = i;
            SETFLOAT (listOut + k, x->x_clusters[i].members[j]);
        };

    selector = gensym ("clusters_list");
    outlet_anything (x->x_listOut, selector, x->x_numInstances, listOut);

    // free memory
    t_freebytes (minDistIdx, 2 * sizeof (t_instanceIdx));

    // free the database memory
    for (i = 0; i < numInstances; i++)
        t_freebytes (clusterData[i].data, clusterData[i].length * sizeof (t_float));

    t_freebytes (clusterData, numInstances * sizeof (t_instance));

    t_freebytes (pairDists, numPairs * sizeof (t_float));
    t_freebytes (listOut, numInstances * sizeof (t_atom));

    post ("%s: instances clustered.", x->x_objSymbol->s_name);

    } // end of main if/else
}


static void timbreID_uncluster (t_timbreID* x)
{
    t_instanceIdx i;

    // free each x->x_clusters list's memory
    for (i = 0; i < x->x_numClusters; i++)
        t_freebytes (x->x_clusters[i].members, x->x_clusters[i].numMembers * sizeof (t_instanceIdx));

// don't need to do this if x_clusters is always x_numInstances in size
// x->x_clusters = (t_cluster *)t_resizebytes (x->x_clusters, x->x_numClusters * sizeof (t_cluster), x->x_numInstances * sizeof (t_cluster));

    for (i = 0; i < x->x_numInstances; i++)
        x->x_clusters[i].members = (t_instanceIdx *)t_getbytes (2 * sizeof (t_instanceIdx));

    x->x_numClusters = x->x_numInstances;

    for (i = 0; i < x->x_numInstances; i++)
    {
        x->x_instances[i].clusterMembership = i; // init membership to index
        x->x_clusters[i].members[0] = i; // first member of the cluster is the instance index
        x->x_clusters[i].members[1] = UINT_MAX; // terminate cluster member list with UINT_MAX
        x->x_clusters[i].numMembers = 2;
        x->x_clusters[i].votes = 0;
    }

    post ("%s: instances unclustered.", x->x_objSymbol->s_name);
}


static void timbreID_computeVariance (t_timbreID* x)
{
    if (x->x_numInstances > 0)
    {
        t_instanceIdx i, j;
        t_float max;
        t_float* attributeVar;
        t_instance* meanCentered;

        // create local memory
        attributeVar = (t_float *)t_getbytes (x->x_maxFeatureLength * sizeof (t_float));
        meanCentered = (t_instance *)t_getbytes (x->x_numInstances * sizeof (t_instance));

        for (i = 0; i < x->x_numInstances; i++)
        {
            meanCentered[i].length = x->x_instances[i].length;
            meanCentered[i].data = (t_float *)t_getbytes (meanCentered[i].length * sizeof (t_float));
        }

        // init mean centered
        for (i = 0; i < x->x_numInstances; i++)
            for (j = 0; j < meanCentered[i].length; j++)
                meanCentered[i].data[j] = 0.0;

        // get the mean of each attribute
        for (i = 0; i < x->x_maxFeatureLength; i++)
            attributeVar[i] = timbreID_attributeMean (x->x_numInstances, i, x->x_instances, x->x_normalize, x->x_attributeData);

        // center the data and write the matrix B
        for (i = 0; i < x->x_numInstances; i++)
            for (j = 0; j < meanCentered[i].length; j++)
            {
                if (x->x_normalize)
                    meanCentered[i].data[j] = ((x->x_instances[i].data[j] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar) - attributeVar[j];
                else
                    meanCentered[i].data[j] = x->x_instances[i].data[j] - attributeVar[j];
            }

        // variance is calculated as: sum (B(:,1).^2)/(M - 1) for the first attribute
        // run process by matrix columns rather than rows, hence the j, i order
        for (j = 0; j < x->x_maxFeatureLength; j++)
        {
            attributeVar[j] = 0;

            for (i = 0; i < x->x_numInstances; i++)
            {
                if (j < meanCentered[i].length)
                    attributeVar[j] += meanCentered[i].data[j] * meanCentered[i].data[j];
            }

            if ((x->x_numInstances - 1) > 0)
                attributeVar[j] /= x->x_numInstances - 1;
        }

        // sort attributeOrder by largest variances: find max in attributeVar,
        // replace it with -FLT_MAX, find next max.
        for (i = 0; i < x->x_maxFeatureLength; i++)
        {
            max = 0.0;

            for (j = 0; j < x->x_maxFeatureLength; j++)
            {
                if (attributeVar[j] > max)
                {
                    max = attributeVar[j];
                    x->x_attributeData[i].order = j;
                }
            };

            attributeVar[x->x_attributeData[i].order] = -FLT_MAX;
        };

        // free local memory
        t_freebytes (attributeVar, x->x_maxFeatureLength * sizeof (t_float));

        // free the meanCentered memory
        for (i = 0; i < x->x_numInstances; i++)
            t_freebytes (meanCentered[i].data, meanCentered[i].length * sizeof (t_float));

        t_freebytes (meanCentered, x->x_numInstances * sizeof (t_instance));

        post ("%s: attributes ordered by variance.", x->x_objSymbol->s_name);
    }
    else
        pd_error (x, "%s: no instances for variance computation.", x->x_objSymbol->s_name);

}


static void timbreID_clustersList (t_timbreID* x)
{
    t_instanceIdx i, j, k;
    t_atom* listOut;
    t_symbol* selector;

    // create local memory
    listOut = (t_atom *)t_getbytes (x->x_numInstances * sizeof (t_atom));

    for (i = 0, k = 0; i < x->x_numClusters; i++)
        for (j = 0; j < x->x_clusters[i].numMembers - 1; j++, k++) // -1 because it's terminated by UINT_MAX
            SETFLOAT (listOut + k, x->x_clusters[i].members[j]);

    selector = gensym ("clusters_list");
    outlet_anything (x->x_listOut, selector, x->x_numInstances, listOut);

    // free local memory
    t_freebytes (listOut, x->x_numInstances * sizeof (t_atom));
}


static void timbreID_clusterList (t_timbreID* x, t_floatarg idx)
{
    t_instanceIdx idxInt;

    idxInt = (idx < 0) ? 0 : idx;

    if (idxInt >= x->x_numClusters)
        pd_error (x, "%s: cluster %i does not exist.", x->x_objSymbol->s_name, idxInt);
    else
    {
        t_instanceIdx i, numMembers;
        t_atom* listOut;
        t_symbol* selector;

        numMembers = x->x_clusters[idxInt].numMembers - 1;

        // create local memory
        listOut = (t_atom *)t_getbytes (numMembers * sizeof (t_atom));

        for (i = 0; i < numMembers; i++)
            SETFLOAT (listOut + i, x->x_clusters[idxInt].members[i]);

        selector = gensym ("cluster_list");
        outlet_anything (x->x_listOut, selector, numMembers, listOut);

        // free local memory
        t_freebytes (listOut, numMembers * sizeof (t_atom));
    }
}


static void timbreID_clusterMembership (t_timbreID* x, t_floatarg idx)
{
    t_instanceIdx idxInt;

    idxInt = (idx < 0) ? 0 : idx;

    if (idxInt >= x->x_numInstances)
        pd_error (x, "%s: instance %i does not exist.", x->x_objSymbol->s_name, idxInt);
    else
    {
        t_atom listOut;
        t_symbol* selector;

        SETFLOAT (&listOut, x->x_instances[idxInt].clusterMembership);
        selector = gensym ("cluster_membership");
        outlet_anything (x->x_listOut, selector, 1, &listOut);
    }
}


static void timbreID_computeOrder (t_timbreID* x, t_floatarg reference)
{
    t_instanceIdx i, j, smallIdx, ref;
    t_float smallest, dist;
    t_instance* instancesCopy;
    t_atom* listOut;
    t_symbol* selector;

    // create local memory
    instancesCopy = (t_instance *)t_getbytes (x->x_numInstances * sizeof (t_instance));
    listOut = (t_atom *)t_getbytes (x->x_numInstances * sizeof (t_atom));

    for (i = 0; i < x->x_numInstances; i++)
        instancesCopy[i].data = (t_float *)t_getbytes (x->x_instances[i].length * sizeof (t_float));

    if (reference > x->x_numInstances - 1)
        ref = x->x_numInstances - 1;
    else if (reference < 0)
        ref = 0;
    else
        ref = reference;

    // make a local copy of instances so they can be abused
    for (i = 0; i < x->x_numInstances; i++)
    {
        instancesCopy[i].length = x->x_instances[i].length;
        instancesCopy[i].clusterMembership = x->x_instances[i].clusterMembership;

        for (j = 0; j < instancesCopy[i].length; j++)
            instancesCopy[i].data[j] = x->x_instances[i].data[j];
    }

    for (i = 0; i < x->x_numInstances; i++)
    {
        smallest = FLT_MAX;
        smallIdx = 0;

        for (j = 0; j < x->x_numInstances; j++)
        {
            // skip this iteration if this instance slot has already been used.
            if (instancesCopy[j].data[0] == FLT_MAX)
                continue;

            dist = timbreID_getDist (x, x->x_instances[ref], instancesCopy[j]);

            if (dist < smallest)
            {
                smallest = dist;
                smallIdx = j;
            };

        };

        SETFLOAT (listOut + i, smallIdx); // store the best from this round;

        if (x->x_relativeOrdering)
            ref = smallIdx; // reorient search to nearest match;

        // set this instance to something huge so it will never be chosen as a good match
        for (j = 0; j < x->x_instances[smallIdx].length; j++)
            instancesCopy[smallIdx].data[j] = FLT_MAX;
    };

    selector = gensym ("order");
    outlet_anything (x->x_listOut, selector, x->x_numInstances, listOut);

    // free local memory
    for (i = 0; i < x->x_numInstances; i++)
        t_freebytes (instancesCopy[i].data, instancesCopy[i].length * sizeof (t_float));

    t_freebytes (instancesCopy, x->x_numInstances * sizeof (t_instance));
    t_freebytes (listOut, x->x_numInstances * sizeof (t_atom));
}


static void timbreID_relativeOrdering (t_timbreID* x, t_floatarg rel)
{
    if (rel < 0)
        x->x_relativeOrdering = false;
    else if (rel > 1)
        x->x_relativeOrdering = true;
    else
        x->x_relativeOrdering = rel;

    if (x->x_relativeOrdering)
        post ("%s: relative ordering ON.", x->x_objSymbol->s_name);
    else
        post ("%s: relative ordering OFF.", x->x_objSymbol->s_name);
}


static void timbreID_distMetric (t_timbreID* x, t_floatarg d)
{
    if (d < 0)
        x->x_distMetric = 0;
    if (d > 3)
        x->x_distMetric = 3;
    else
        x->x_distMetric = d;

    switch (x->x_distMetric)
    {
        case euclidean:
            post ("%s: distance metric: EUCLIDEAN.", x->x_objSymbol->s_name);
            break;
        case taxi:
            post ("%s: distance metric: MANHATTAN (taxicab distance).", x->x_objSymbol->s_name);
            break;
        case correlation:
            post ("%s: distance metric: PEARSON CORRELATION COEFF.", x->x_objSymbol->s_name);
            break;
        default:
            break;
    };
}


static void timbreID_weights (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    t_attributeIdx i, listLength;

    listLength = argc;

    if (listLength > x->x_maxFeatureLength)
    {
        pd_error (x, "%s: weights list longer than current feature length. ignoring excess weights", x->x_objSymbol->s_name);
        argc = x->x_maxFeatureLength;
    }

    for (i = 0; i < listLength; i++)
        x->x_attributeData[i].weight = atom_getfloat (argv + i);

    // if only the first few of a long feature vector are specified, fill in the rest with 1.0
    for (; i < x->x_maxFeatureLength; i++)
        x->x_attributeData[i].weight = 1.0;
}


static void timbreID_attributeNames (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    t_attributeIdx i, listLength;

    listLength = argc;

    if (listLength > x->x_maxFeatureLength)
    {
        pd_error (x, "%s: attribute list longer than current max feature length. ignoring excess attributes", x->x_objSymbol->s_name);
        listLength = x->x_maxFeatureLength;
    }

    for (i = 0; i < listLength; i++)
        x->x_attributeData[i].name = atom_getsymbol (argv + i);

    post ("%s: attribute names received.", x->x_objSymbol->s_name);
}


static void timbreID_attributeInfo (t_timbreID* x, t_floatarg idx)
{
    t_attributeIdx thisIdx;
    thisIdx = (idx < 0) ? 0 : idx;

    if (thisIdx >= x->x_maxFeatureLength)
        pd_error (x, "%s: attribute %i does not exist", x->x_objSymbol->s_name, thisIdx);
    else
    {
        t_symbol* selector;
        t_atom listOut[3];

        selector = gensym ("attribute_info");
        // name, weight, order
        SETSYMBOL(&listOut[0], x->x_attributeData[thisIdx].name);
        SETFLOAT (&listOut[1], x->x_attributeData[thisIdx].weight);
        SETFLOAT (&listOut[2], x->x_attributeData[thisIdx].order);
        outlet_anything (x->x_listOut, selector, 3, &listOut[0]);
    }
}


static void timbreID_attributeOrder (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    t_attributeIdx i, j, listLength;

    listLength = argc;

    if (listLength > x->x_maxFeatureLength)
    {
        pd_error (x, "%s: attribute list longer than current max feature length. ignoring excess attributes", x->x_objSymbol->s_name);
        listLength = x->x_maxFeatureLength;
    }

    for (i = 0; i < listLength; i++)
    {
        x->x_attributeData[i].order = atom_getfloat (argv + i);

        if (x->x_attributeData[i].order > x->x_maxFeatureLength - 1)
        {
            // initialize attributeOrder
            for (j = 0; j < x->x_maxFeatureLength; j++)
                x->x_attributeData[j].order = j;

            pd_error (x, "%s: attribute %i out of range. attribute order initialized.", x->x_objSymbol->s_name, x->x_attributeData[i].order);
        }
    }

    // fill any remainder with attribute 0
    for (; i < x->x_maxFeatureLength; i++)
        x->x_attributeData[i].order = 0;

    post ("%s: attributes re-ordered.", x->x_objSymbol->s_name);
}


static void timbreID_attributeRange (t_timbreID* x, t_floatarg lo, t_floatarg hi)
{
    x->x_attributeLo = (lo < 0) ? 0 : lo;
    x->x_attributeLo = (x->x_attributeLo >= x->x_maxFeatureLength) ? x->x_maxFeatureLength - 1 : x->x_attributeLo;

    x->x_attributeHi = (hi < 0) ? 0 : hi;
    x->x_attributeHi = (x->x_attributeHi >= x->x_maxFeatureLength) ? x->x_maxFeatureLength - 1 : x->x_attributeHi;

    if (x->x_attributeLo > x->x_attributeHi)
    {
        t_attributeIdx tmp;
        tmp = x->x_attributeHi;
        x->x_attributeHi = x->x_attributeLo;
        x->x_attributeLo = tmp;
    }

    post ("%s: attribute range: %i through %i.", x->x_objSymbol->s_name, x->x_attributeLo, x->x_attributeHi);
}


static void timbreID_reorderAttributes (t_timbreID* x)
{
    t_attributeIdx i;

    // initialize attributeOrder
    for (i = 0; i < x->x_maxFeatureLength; i++)
        x->x_attributeData[i].order = i;

    post ("%s: attribute order initialized.", x->x_objSymbol->s_name);
}


static void timbreID_featureList (t_timbreID* x, t_floatarg idx, t_floatarg normFlag)
{
    t_instanceIdx idxInt;
    t_bool normalize, normRange;

    post ("%s WARNING: \"feature_list\" method name is deprecated, use \"instance_list\" instead", x->x_objSymbol->s_name);

    normalize = false;
    normRange = false;

    switch ((int)normFlag)
    {
        case 0:
            normalize = false;
            normRange = false;
            break;
        case 1:
            normalize = true;
            normRange = false;
            break;
        case 2:
            normalize = true;
            normRange = true;
            break;
        default:
            break;
    }

    idxInt = (idx < 0) ? 0 : idx;

    if (idxInt >= x->x_numInstances)
        pd_error (x, "%s: instance %i does not exist.", x->x_objSymbol->s_name, idxInt);
    else
    {
        t_attributeIdx i, thisFeatureLength;
        t_atom* listOut;
        t_symbol* selector;

        thisFeatureLength = x->x_instances[idxInt].length;

        normalize = (normalize < 0) ? 0 : normalize;
        normalize = (normalize > 1) ? 1 : normalize;

        normRange = (normRange < 0) ? 0 : normRange;
        normRange = (normRange > 1) ? 1 : normRange;

        // create local memory
        listOut = (t_atom *)t_getbytes (thisFeatureLength * sizeof (t_atom));

        for (i = 0; i < thisFeatureLength; i++)
        {
            if (normalize)
            {
                if ( !x->x_normalize)
                {
                    pd_error (x, "%s: feature database not normalized. cannot output normalized values.", x->x_objSymbol->s_name);

                    // free local memory before exit
                    t_freebytes (listOut, thisFeatureLength * sizeof (t_atom));

                    return;
                }

                if (normRange)
                    SETFLOAT (listOut + i, ((x->x_instances[idxInt].data[i] - x->x_attributeData[i].normData.minVal) * x->x_attributeData[i].normData.normScalar * 2.0) - 1.0);
                else
                    SETFLOAT (listOut + i, (x->x_instances[idxInt].data[i] - x->x_attributeData[i].normData.minVal) * x->x_attributeData[i].normData.normScalar);
            }
            else
                SETFLOAT (listOut + i, x->x_instances[idxInt].data[i]);
        }

        selector = gensym ("feature_list");
        outlet_anything (x->x_listOut, selector, thisFeatureLength, listOut);

        // free local memory
        t_freebytes (listOut, thisFeatureLength * sizeof (t_atom));
    }
}


static void timbreID_instanceList (t_timbreID* x, t_floatarg idx, t_floatarg normFlag)
{
    t_instanceIdx idxInt;
    t_bool normalize, normRange;

    normalize = false;
    normRange = false;

    switch ((int)normFlag)
    {
        case 0:
            normalize = false;
            normRange = false;
            break;
        case 1:
            normalize = true;
            normRange = false;
            break;
        case 2:
            normalize = true;
            normRange = true;
            break;
        default:
            break;
    }

    idxInt = (idx < 0) ? 0 : idx;

    if (idxInt >= x->x_numInstances)
        pd_error (x, "%s: instance %i does not exist.", x->x_objSymbol->s_name, idxInt);
    else
    {
        t_attributeIdx i, thisFeatureLength;
        t_atom* listOut;
        t_symbol* selector;

        thisFeatureLength = x->x_instances[idxInt].length;

        normalize = (normalize < 0) ? 0 : normalize;
        normalize = (normalize > 1) ? 1 : normalize;

        normRange = (normRange < 0) ? 0 : normRange;
        normRange = (normRange > 1) ? 1 : normRange;

        // create local memory
        listOut = (t_atom *)t_getbytes (thisFeatureLength * sizeof (t_atom));

        for (i = 0; i < thisFeatureLength; i++)
        {
            if (normalize)
            {
                if ( !x->x_normalize)
                {
                    pd_error (x, "%s: feature database not normalized. cannot output normalized values.", x->x_objSymbol->s_name);

                    // free local memory before exit
                    t_freebytes (listOut, thisFeatureLength * sizeof (t_atom));

                    return;
                }

                if (normRange)
                    SETFLOAT (listOut + i, ((x->x_instances[idxInt].data[i] - x->x_attributeData[i].normData.minVal) * x->x_attributeData[i].normData.normScalar * 2.0) - 1.0);
                else
                    SETFLOAT (listOut + i, (x->x_instances[idxInt].data[i] - x->x_attributeData[i].normData.minVal) * x->x_attributeData[i].normData.normScalar);
            }
            else
                SETFLOAT (listOut + i, x->x_instances[idxInt].data[i]);
        }

        selector = gensym ("instance_list");
        outlet_anything (x->x_listOut, selector, thisFeatureLength, listOut);

        // free local memory
        t_freebytes (listOut, thisFeatureLength * sizeof (t_atom));
    }
}


static void timbreID_attributeList (t_timbreID* x, t_floatarg idx, t_floatarg normFlag)
{
    t_attributeIdx idxInt;
    t_bool normalize, normRange;

    normalize = false;
    normRange = false;

    switch ((int)normFlag)
    {
        case 0:
            normalize = false;
            normRange = false;
            break;
        case 1:
            normalize = true;
            normRange = false;
            break;
        case 2:
            normalize = true;
            normRange = true;
            break;
        default:
            break;
    }

    idxInt = (idx < 0) ? 0 : idx;

    if (x->x_numInstances == 0 || idxInt > x->x_maxFeatureLength - 1)
        pd_error (x, "%s: attribute %i does not exist.", x->x_objSymbol->s_name, idxInt);
    else
    {
        t_attributeIdx i, attributeListLength;
        t_atom* listOut;
        t_symbol* selector;

        attributeListLength = x->x_numInstances;

        normalize = (normalize < 0) ? 0 : normalize;
        normalize = (normalize > 1) ? 1 : normalize;

        normRange = (normRange < 0) ? 0 : normRange;
        normRange = (normRange > 1) ? 1 : normRange;

        // create local memory
        listOut = (t_atom *)t_getbytes (attributeListLength * sizeof (t_atom));

        for (i = 0; i < attributeListLength; i++)
        {
            if (idxInt >= x->x_instances[i].length)
            {
                SETFLOAT (listOut + i, 0.0);
                post ("%s WARNING: attribute %i does not exist for instance %i, outputting 0 instead.", x->x_objSymbol->s_name, idxInt, i);
            }
            else
            {
                if (normalize)
                {
                    if ( !x->x_normalize)
                    {
                        pd_error (x, "%s: feature database not normalized. cannot output normalized values.", x->x_objSymbol->s_name);

                        // free local memory before exit
                        t_freebytes (listOut, attributeListLength * sizeof (t_atom));

                        return;
                    }

                    if (normRange)
                        SETFLOAT (listOut + i, ((x->x_instances[i].data[idxInt] - x->x_attributeData[idxInt].normData.minVal) * x->x_attributeData[idxInt].normData.normScalar * 2.0) - 1.0);
                    else
                        SETFLOAT (listOut + i, (x->x_instances[i].data[idxInt] - x->x_attributeData[idxInt].normData.minVal) * x->x_attributeData[idxInt].normData.normScalar);
                }
                else
                    SETFLOAT (listOut + i, x->x_instances[i].data[idxInt]);
            }
        }

        selector = gensym ("attribute_list");
        outlet_anything (x->x_listOut, selector, attributeListLength, listOut);

        // free local memory
        t_freebytes (listOut, attributeListLength * sizeof (t_atom));
    }
}


static void timbreID_similarityMatrix (t_timbreID* x, t_floatarg startInstance, t_floatarg finishInstance, t_floatarg normalize)
{
    if (x->x_numInstances)
    {
        t_instanceIdx i, j, k, l, numInst, startInst, finishInst;
        t_bool norm;

        startInstance = (startInstance < 0) ? 0 : startInstance;
        startInstance = (startInstance >= x->x_numInstances) ? x->x_numInstances - 1 : startInstance;
        startInst = startInstance;

        finishInstance = (finishInstance < 0) ? 0 : finishInstance;
        finishInstance = (finishInstance >= x->x_numInstances) ? x->x_numInstances - 1 : finishInstance;
        finishInst = finishInstance;

        normalize = (normalize < 0) ? 0 : normalize;
        normalize = (normalize > 1) ? 1 : normalize;
        norm = normalize;

        if (startInst > finishInst)
        {
            t_instanceIdx tmp;

            tmp = finishInst;
            finishInst = startInst;
            startInst = tmp;
        }

        if (finishInst - startInst == 0)
        {
            pd_error (x, "%s: bad instance range for similarity matrix", x->x_objSymbol->s_name);
            return;
        }

        numInst = finishInst - startInst + 1;

        if (numInst)
        {
            t_float maxDist;
            t_instance* distances;
            t_symbol* selector;

            // create local memory
            // only using t_instance type as a convenient way to make a numInst X numInst matrix
            distances = (t_instance *)t_getbytes (numInst * sizeof (t_instance));

            for (i = 0; i < numInst; i++)
            {
                distances[i].data = (t_float *)t_getbytes (numInst * sizeof (t_float));
                distances[i].length = numInst;

                for (j = 0; j < distances[i].length; j++)
                    distances[i].data[j] = 0.0;
            }

            maxDist = -1;

            for (i = startInst, j = 0; i <= finishInst; i++, j++)
            {
                for (k = startInst, l = 0; k <= finishInst; k++, l++)
                {
                    t_float dist;

                    dist = timbreID_getDist (x, x->x_instances[i], x->x_instances[k]);

                    if (dist > maxDist)
                        maxDist = dist;

                    if (l < distances[i].length)
                        distances[j].data[l] = dist;
                }
            }

            maxDist = 1.0 / maxDist;

            for (i = startInst; i < numInst; i++)
            {
                t_atom* listOut;

                listOut = (t_atom *)t_getbytes (numInst * sizeof (t_atom));

                for (j = 0; j < distances[i].length; j++)
                {
                    t_float dist;

                    dist = distances[i].data[j];

                    if (norm)
                        dist *= maxDist;

                    SETFLOAT (listOut + j, dist);
                }

                selector = gensym ("similarity_matrix");
                outlet_anything (x->x_listOut, selector, numInst, listOut);

                t_freebytes (listOut, numInst * sizeof (t_atom));
            }


            // free local memory
            for (i = 0; i < numInst; i++)
                t_freebytes (distances[i].data, numInst * sizeof (t_float));

            t_freebytes (distances, numInst * sizeof (t_instance));

        }
        else
            pd_error (x, "%s: bad range of instances", x->x_objSymbol->s_name);
    }
    else
        pd_error (x, "%s: no training instances have been loaded.", x->x_objSymbol->s_name);
}


static void timbreID_minValues (t_timbreID* x)
{
    if (x->x_normalize)
    {
        t_attributeIdx i;
        t_atom* listOut;
        t_symbol* selector;

        // create local memory
        listOut = (t_atom *)t_getbytes (x->x_maxFeatureLength * sizeof (t_atom));

        for (i = 0; i < x->x_maxFeatureLength; i++)
            SETFLOAT (listOut + i, x->x_attributeData[i].normData.minVal);

        selector = gensym ("min_values");
        outlet_anything (x->x_listOut, selector, x->x_maxFeatureLength, listOut);

        // free local memory
        t_freebytes (listOut, x->x_maxFeatureLength * sizeof (t_atom));
    }
    else
        pd_error (x, "%s: feature database not normalized. minimum values not calculated yet.", x->x_objSymbol->s_name);
}


static void timbreID_maxValues (t_timbreID* x)
{
    if (x->x_normalize)
    {
        t_attributeIdx i;
        t_atom* listOut;
        t_symbol* selector;

        // create local memory
        listOut = (t_atom *)t_getbytes (x->x_maxFeatureLength * sizeof (t_atom));

        for (i = 0; i < x->x_maxFeatureLength; i++)
            SETFLOAT (listOut + i, x->x_attributeData[i].normData.maxVal);

        selector = gensym ("max_values");
        outlet_anything (x->x_listOut, selector, x->x_maxFeatureLength, listOut);

        // free local memory
        t_freebytes (listOut, x->x_maxFeatureLength * sizeof (t_atom));
    }
    else
        pd_error (x, "%s: feature database not normalized. maximum values not calculated yet.", x->x_objSymbol->s_name);
}


static void timbreID_mink (t_timbreID* x, t_floatarg k, t_floatarg normFlag)
{
    t_sampIdx i;
    t_attributeIdx j, kInt;
    t_float* attVals;
    t_symbol* selector;
    t_atom* outputList;
    t_bool normalize, normRange;

    if (k < 0 || k >= x->x_numInstances)
    {
        pd_error (x, "%s: k must be less than the number of instances in the database.", x->x_objSymbol->s_name);
        // TODO: shouldn't have an else branch after a return statement
        return;
    }
    else
        kInt = k;

    normalize = false;
    normRange = false;

    switch ((int)normFlag)
    {
        case 0:
            normalize = false;
            normRange = false;
            break;
        case 1:
            normalize = true;
            normRange = false;
            break;
        case 2:
            normalize = true;
            normRange = true;
            break;
        default:
            break;
    }

    attVals = (t_float *)t_getbytes (x->x_numInstances * sizeof (t_float));
    outputList = (t_atom *)t_getbytes ((kInt + 1) * sizeof (t_atom));

    selector = gensym ("min_k_values");

    for (j = 0; j < x->x_maxFeatureLength; j++)
    {
        for (i = 0; i < x->x_numInstances; i++)
            attVals[i] = x->x_instances[i].data[j];

        tIDLib_bubbleSort (x->x_numInstances, attVals);

        // start the list with the ID of this attribute
        SETFLOAT (outputList, j);

        for (i = 0; i < kInt; i++)
        {
            if (normalize)
            {
                if ( !x->x_normalize)
                {
                    pd_error (x, "%s: feature database not normalized. cannot output normalized values.", x->x_objSymbol->s_name);

                    // free local memory before exit
                    t_freebytes (attVals, x->x_numInstances * sizeof (t_float));
                    t_freebytes (outputList, (kInt + 1) * sizeof (t_atom));

                    return;
                }

                if (normRange)
                    SETFLOAT (outputList + 1 + i, ((attVals[i] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar * 2.0) - 1.0);
                else
                    SETFLOAT (outputList + 1 + i, (attVals[i] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar);
            }
            else
                SETFLOAT (outputList + 1 + i, attVals[i]);
        }

        // output k values for this attribute
        outlet_anything (x->x_listOut, selector, kInt + 1, outputList);
    }

    // free the attVals and outputList buffers
    t_freebytes (attVals, x->x_numInstances * sizeof (t_float));
    t_freebytes (outputList, (kInt + 1) * sizeof (t_atom));
}


static void timbreID_maxk (t_timbreID* x, t_floatarg k, t_floatarg normFlag)
{
    t_sampIdx i;
    t_attributeIdx j, kInt;
    t_float* attVals;
    t_symbol* selector;
    t_atom* outputList;
    t_bool normalize, normRange;

    if (k < 0 || k >= x->x_numInstances)
    {
        pd_error (x, "%s: k must be less than the number of instances in the database.", x->x_objSymbol->s_name);
        // TODO: shouldn't have an else branch after a return statement
        return;
    }
    else
        kInt = k;

    normalize = false;
    normRange = false;

    switch ((int)normFlag)
    {
        case 0:
            normalize = false;
            normRange = false;
            break;
        case 1:
            normalize = true;
            normRange = false;
            break;
        case 2:
            normalize = true;
            normRange = true;
            break;
        default:
            break;
    }

    attVals = (t_float *)t_getbytes (x->x_numInstances * sizeof (t_float));
    outputList = (t_atom *)t_getbytes ((kInt + 1) * sizeof (t_atom));

    selector = gensym ("max_k_values");

    for (j = 0; j < x->x_maxFeatureLength; j++)
    {
        for (i = 0; i < x->x_numInstances; i++)
            attVals[i] = x->x_instances[i].data[j];

        tIDLib_bubbleSort (x->x_numInstances, attVals);

        // start the list with the ID of this attribute
        SETFLOAT (outputList, j);

        for (i = 0; i < kInt; i++)
        {
            if (normalize)
            {
                if ( !x->x_normalize)
                {
                    pd_error (x, "%s: feature database not normalized. cannot output normalized values.", x->x_objSymbol->s_name);

                    // free local memory before exit
                    t_freebytes (attVals, x->x_numInstances * sizeof (t_float));
                    t_freebytes (outputList, (kInt + 1) * sizeof (t_atom));

                    return;
                }

                if (normRange)
                    SETFLOAT (outputList + 1 + i, ((attVals[(x->x_numInstances - 1) - i] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar * 2.0) - 1.0);
                else
                    SETFLOAT (outputList + 1 + i, (attVals[(x->x_numInstances - 1) - i] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar);
            }
            else
                SETFLOAT (outputList + 1 + i, attVals[(x->x_numInstances - 1) - i]);
        }

        // output k values for this attribute
        outlet_anything (x->x_listOut, selector, kInt + 1, outputList);
    }

    // free the attVals and outputList buffers
    t_freebytes (attVals, x->x_numInstances * sizeof (t_float));
    t_freebytes (outputList, (kInt + 1) * sizeof (t_atom));
}


static void timbreID_classReference (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    t_symbol* commandName;

    if (x->x_numInstances == 0)
    {
        pd_error (x, "%s: no training instances have been loaded. cannot define a class reference.", x->x_objSymbol->s_name);
        return;
    }

    commandName = atom_getsymbol (argv);

    if ( !strcmp (commandName->s_name, "define"))
    {
        t_instanceIdx classIdx, i;
        t_attributeIdx numAttributes, j;
        t_float* attVals;

        // within this "define" command, could offer the option of using min/max K values approach vs making average vectors based on all instances in each cluster (would require clustering in advance).

        // store the class index
        classIdx = x->x_numClassRefs;

        // the number of attributes in this definition is the number of remaining arguments after the "define" symbol
        numAttributes = argc - 1;
        // post ("numAttributes: %i", numAttributes);

        // confirm that the definition specifies a K percent for all attributes
        // if not, return with an error before we do any memory management below
        if (numAttributes != x->x_minFeatureLength)
        {
            pd_error (x, "%s: class reference definition must provide a K percent value for all attributes. failed to define class reference.", x->x_objSymbol->s_name);
            return;
        }

        // add an instance to the classRefs database
        x->x_classRefs = (t_instance *)t_resizebytes (x->x_classRefs, x->x_numClassRefs * sizeof (t_instance), (x->x_numClassRefs + 1) * sizeof (t_instance));

        // update the number of class references
        x->x_numClassRefs++;

        // store the number of attributes in the .length field
        x->x_classRefs[classIdx].length = numAttributes;

        // get memory for the actual data
        x->x_classRefs[classIdx].data = (t_float *)t_getbytes (numAttributes * sizeof (t_float));

        // get local memory to store each attribute column
        attVals = (t_float *)t_getbytes (x->x_numInstances * sizeof (t_float));

        // go through each attribute and:
            // get the average of the min/max K values of the specified attribute
        for (j = 0; j < numAttributes; j++)
        {
            t_float pct, kSum;
            t_instanceIdx kInstances;

            pct = atom_getfloat (argv + 1 + j);

            // get this attribute column
            for (i = 0; i < x->x_numInstances; i++)
                attVals[i] = x->x_instances[i].data[j];

            // sort it
            tIDLib_bubbleSort (x->x_numInstances, attVals);

            // convert pct into a rounded K value
            kInstances = (fabs (pct) / 100.0) * x->x_numInstances + 0.5;
            // post ("kInstances[%i]: %i", j, kInstances);

            // initialize result, where we'll store the average
            kSum = 0.0;

            // if pct is negative, run mink
            // if pct is positive, run maxk
            if (pct < 0)
            {
                for (i = 0; i < kInstances; i++)
                {
                    if (x->x_normalize)
                        kSum += (attVals[i] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar;
                    else
                        kSum += attVals[i];
                }
            }
            else
            {
                for (i = 0; i < kInstances; i++)
                {
                    if (x->x_normalize)
                        kSum += (attVals[(x->x_numInstances - 1) - i] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar;
                    else
                        kSum += attVals[(x->x_numInstances - 1) - i];
                }
            }

            // divide the running sum in kSum by K to get the average of the min K values
            x->x_classRefs[classIdx].data[j] = kSum / kInstances;
            // post ("avg[%i]: %f", j, x->x_classRefs[classIdx].data[j]);
        }

        // free the attVals buffer
        t_freebytes (attVals, x->x_numInstances * sizeof (t_float));

        post ("%s: defined class reference %i.", x->x_objSymbol->s_name, classIdx);
    }
    else if ( !strcmp (commandName->s_name, "clear"))
    {
        t_instanceIdx i;

        // free the class reference memory per instance
        for (i = 0; i < x->x_numClassRefs; i++)
            t_freebytes (x->x_classRefs[i].data, x->x_classRefs[i].length * sizeof (t_float));

        // resize the class ref database to 0 bytes
        x->x_classRefs = (t_instance *)t_resizebytes (x->x_classRefs, x->x_numClassRefs * sizeof (t_instance), 0);

        x->x_numClassRefs = 0;

        post ("%s: cleared all class reference vectors.", x->x_objSymbol->s_name);
    }
    else if ( !strcmp (commandName->s_name, "id"))
    {
        t_attributeIdx vecLen;
        t_instanceIdx i, winningIdx;
        t_float* inputVector;
        t_float winningDist;

        if (x->x_numClassRefs == 0)
        {
            pd_error (x, "%s: no class references have been defined.", x->x_objSymbol->s_name);
            return;
        }

        //determine the length of the input vector (one less than argc to remove the "id" command)
        vecLen = argc - 1;

        // get memory for the input vector
        inputVector = (t_float *)t_getbytes (vecLen * sizeof (t_float));

        // fill the inputVector memory with the input data
        // we can't normalize the input vector here since we don't know which attributes are in the class reference. timbreID help file must state clearly that incoming "class_reference id" vectors must be normalized ahead of time
        for (i = 0; i < vecLen; i++)
            inputVector[i] = atom_getfloat (argv + 1 + i);

        winningDist = FLT_MAX;
        winningIdx = 0;

        // compare input vector with all defined class references
        for (i = 0; i < x->x_numClassRefs; i++)
        {
            t_float thisDistance;

            thisDistance = timbreID_getClassRefDist (x, vecLen, inputVector, x->x_classRefs[i]);

            if (thisDistance < winningDist)
            {
                winningDist = thisDistance;
                winningIdx = i;
            }
        }

        // output the winning distance
        outlet_float (x->x_nearestDist, winningDist);

        // output the winning index
        outlet_float (x->x_id, winningIdx);

        // free the inputVector buffer
        t_freebytes (inputVector, vecLen * sizeof (t_float));
    }
    else if ( !strcmp (commandName->s_name, "get"))
    {
        t_instanceIdx i;
        t_attributeIdx j, classVecLen;
        t_symbol* selector;
        t_atom* listOut;

        classVecLen = x->x_classRefs[0].length;

        // make local memory to output the list
        listOut = t_getbytes ((classVecLen + 1) * sizeof (t_atom));

        // fill the list with data from each class reference and output
        for (i = 0; i < x->x_numClassRefs; i++)
        {
            SETFLOAT (listOut, i);

            for (j = 0; j < classVecLen; j++)
                SETFLOAT (listOut + 1 + j, x->x_classRefs[i].data[j]);

            selector = gensym ("class_reference");
            outlet_anything (x->x_listOut, selector, classVecLen + 1, listOut);
        }

        // free local memory
        t_freebytes (listOut, (classVecLen + 1) * sizeof (t_atom));
    }
    else
    {
        pd_error (x, "%s: no such command (%s) for class_reference.", x->x_objSymbol->s_name, commandName->s_name);
    }
}


static void timbreID_minFeatureLength (t_timbreID* x)
{
    if (x->x_numInstances)
    {
        t_atom listOut;
        t_symbol* selector;

        SETFLOAT (&listOut, x->x_minFeatureLength);

        selector = gensym ("min_feature_length");
        outlet_anything (x->x_listOut, selector, 1, &listOut);
    }
    else
        pd_error (x, "%s: no training instances have been loaded.", x->x_objSymbol->s_name);
}


static void timbreID_maxFeatureLength (t_timbreID* x)
{
    if (x->x_numInstances)
    {
        t_atom listOut;
        t_symbol* selector;

        SETFLOAT (&listOut, x->x_maxFeatureLength);

        selector = gensym ("max_feature_length");
        outlet_anything (x->x_listOut, selector, 1, &listOut);
    }
    else
        pd_error (x, "%s: no training instances have been loaded.", x->x_objSymbol->s_name);
}


static void timbreID_clear (t_timbreID* x)
{
    t_instanceIdx i;

    // uncluster first
    timbreID_uncluster (x);

    // free the each instance's data memory
    for (i = 0; i < x->x_numInstances; i++)
        t_freebytes (x->x_instances[i].data, x->x_instances[i].length * sizeof (t_float));

    x->x_instances = (t_instance *)t_resizebytes (x->x_instances, x->x_numInstances * sizeof (t_instance), 0);

    // free each cluster's members memory
    for (i = 0; i < x->x_numClusters; i++)
        t_freebytes (x->x_clusters[i].members, x->x_clusters[i].numMembers * sizeof (t_instanceIdx));

    x->x_clusters = (t_cluster *)t_resizebytes (x->x_clusters, x->x_numInstances * sizeof (t_cluster), 0);

    // postFlag argument FALSE here
    timbreID_attributeDataResize (x, x->x_maxFeatureLength, 0, false);

    x->x_numInstances = 0;
    x->x_numClusters = 0;
    x->x_neighborhood = 0;
    x->x_searchCenter = 0;

    x->x_maxFeatureLength = 0;
    x->x_minFeatureLength = INT_MAX;
    x->x_normalize = false;

    post ("%s: all instances cleared.", x->x_objSymbol->s_name);
}


static void timbreID_write (t_timbreID* x, t_symbol* s)
{
    FILE* filePtr;
    t_instanceIdx i;
    t_instanceIdx* header;
    char fileNameBuf[MAXPDSTRING];

    // make a buffer for the header data, which is the number of instances followed by the length of each instance
    header = (t_instanceIdx *)t_getbytes ((x->x_numInstances + 1) * sizeof (t_instanceIdx));

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "wb");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        // free header memory before exit
        t_freebytes (header, (x->x_numInstances + 1) * sizeof (t_instanceIdx));
        return;
    }

    // write the header information, which is the number of instances, and then the length of each instance
    header[0] = x->x_numInstances;

    for (i = 0; i < x->x_numInstances; i++)
        header[i + 1] = x->x_instances[i].length;

    fwrite (header, sizeof (t_instanceIdx), x->x_numInstances + 1, filePtr);

    // write the instance data
    for (i = 0; i < x->x_numInstances; i++)
        fwrite (x->x_instances[i].data, sizeof (t_float), x->x_instances[i].length, filePtr);

    post ("%s: wrote %i non-normalized instances to %s.", x->x_objSymbol->s_name, x->x_numInstances, fileNameBuf);

    fclose (filePtr);

    // free header memory
    t_freebytes (header, (x->x_numInstances + 1) * sizeof (t_instanceIdx));
}


static void timbreID_read (t_timbreID* x, t_symbol* s)
{
    FILE* filePtr;
    t_instanceIdx i, maxLength, minLength;
    char fileNameBuf[MAXPDSTRING];
    int cTest;
    t_bool txtFlag;

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "rb");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to open %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    // test to see if it's a .timid or .txt file. if any fgetc() results are > 127, it's not a text file.
    txtFlag = true;

    while ((cTest = fgetc (filePtr)) != EOF)
    {
        if (cTest > 127)
            txtFlag = false;
    }

    if (txtFlag)
    {
        pd_error (x, "%s: data in file %s is plain text (.txt), not binary (.timid). use read_text instead.", x->x_objSymbol->s_name, fileNameBuf);
        fclose (filePtr);
        // TODO: shouldn't have an else branch after a return statement
        return;
    }
    else
    {
        // if we're safe, close and re-open so we can get the binary data.
        fclose (filePtr);
        filePtr = fopen (fileNameBuf, "rb");
    }

    maxLength = 0;
    minLength = INT_MAX;

    // erase old instances & clusters and resize to 0. this also does a sub-call to timbreID_attributeDataResize()
    timbreID_clear (x);

    // first item in the header is the number of instances
    fread (&x->x_numInstances, sizeof (t_instanceIdx), 1, filePtr);

    // resize instances & clusterMembers to numInstances
    x->x_instances = (t_instance *)t_resizebytes (x->x_instances, 0, x->x_numInstances * sizeof (t_instance));

    x->x_clusters = (t_cluster *)t_resizebytes (x->x_clusters, 0, x->x_numInstances * sizeof (t_cluster));

    for (i = 0; i < x->x_numInstances; i++)
    {
        // get the length of each instance
        fread (&x->x_instances[i].length, sizeof (t_attributeIdx), 1, filePtr);

        if (x->x_instances[i].length > maxLength)
            maxLength = x->x_instances[i].length;

        if (x->x_instances[i].length < minLength)
            minLength = x->x_instances[i].length;

        // get the appropriate number of bytes for the data
        x->x_instances[i].data = (t_float *)t_getbytes (x->x_instances[i].length * sizeof (t_float));
    }

    x->x_minFeatureLength = minLength;
    x->x_maxFeatureLength = maxLength;
    x->x_neighborhood = x->x_numInstances;
    x->x_numClusters = x->x_numInstances;

    // update x_attributeData based on new x_maxFeatureLength. turn postFlag argument TRUE
    timbreID_attributeDataResize (x, 0, x->x_maxFeatureLength, 1);

    // after loading a database, instances are unclustered
    for (i = 0; i < x->x_numInstances; i++)
    {
        x->x_clusters[i].numMembers = 2;
        x->x_clusters[i].members = (t_instanceIdx *)t_getbytes (x->x_clusters[i].numMembers * sizeof (t_instanceIdx));

        x->x_clusters[i].members[0] = i; // first member of the cluster is the instance index
        x->x_clusters[i].members[1] = UINT_MAX; // terminate with UINT_MAX

        x->x_instances[i].clusterMembership = i; // init instance's cluster membership to index
    }

    // finally, read in the instance data
    for (i = 0; i < x->x_numInstances; i++)
        fread (x->x_instances[i].data, sizeof (t_float), x->x_instances[i].length, filePtr);

    post ("%s: read %i instances from %s.", x->x_objSymbol->s_name, x->x_numInstances, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_writeText (t_timbreID* x, t_symbol* s, t_floatarg normRange)
{
    FILE* filePtr;
    t_instanceIdx i;
    t_attributeIdx j;
    t_float* featurePtr;
    char fileNameBuf[MAXPDSTRING];

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "w");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    normRange = (normRange < 0) ? 0 : normRange;
    normRange = (normRange > 1) ? 1 : normRange;

    for (i = 0; i < x->x_numInstances; i++)
    {
        featurePtr = x->x_instances[i].data;

        j = x->x_attributeLo;

        while (j <= x->x_attributeHi)
        {
            t_float thisFeatureData;
            t_attributeIdx thisAttribute;

            thisAttribute = x->x_attributeData[j].order;

            if (thisAttribute >= x->x_instances[i].length)
                break;

            if (x->x_normalize)
            {
                if (normRange)
                    thisFeatureData = ((*(featurePtr + thisAttribute) - x->x_attributeData[thisAttribute].normData.minVal) * x->x_attributeData[thisAttribute].normData.normScalar * 2.0) - 1;
                else
                    thisFeatureData = (*(featurePtr + thisAttribute) - x->x_attributeData[thisAttribute].normData.minVal) * x->x_attributeData[thisAttribute].normData.normScalar;
            }
            else
                thisFeatureData = *(featurePtr + thisAttribute);

            // What's the best float resolution to print?
            // no space after final value on an instance line
            if (j == x->x_attributeHi)
                fprintf (filePtr, "%0.6f", thisFeatureData);
            else
                fprintf (filePtr, "%0.6f ", thisFeatureData);

            j++;
        };

        fprintf (filePtr, "\n");
       };

    post ("%s: wrote %i instances to %s.", x->x_objSymbol->s_name, x->x_numInstances, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_readText (t_timbreID* x, t_symbol* s)
{
    FILE* filePtr;
    t_instanceIdx i, j, numInstances, stringLength, numSpaces, minLength, maxLength;
    t_float* featurePtr;
    char fileNameBuf[MAXPDSTRING], textLine[TID_MAXTEXTSTRING];

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "r");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to open %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    numInstances = 0;

    while (fgets (textLine, TID_MAXTEXTSTRING, filePtr))
        numInstances++;

    maxLength = 0;
    minLength = INT_MAX;

    // erase old instances & clusters and resize to 0. this also does a sub-call to timbreID_attributeDataResize()
    timbreID_clear (x);

    x->x_numInstances = numInstances;

    // resize instances & clusterMembers to numInstances
    x->x_instances = (t_instance *)t_resizebytes (x->x_instances, 0, x->x_numInstances * sizeof (t_instance));

    x->x_clusters = (t_cluster *)t_resizebytes (x->x_clusters, 0, x->x_numInstances * sizeof (t_cluster));

    // now that we have numInstances, close and re-open file to get the length of each line
    fclose (filePtr);
    filePtr = fopen (fileNameBuf, "r");

    i = 0;

    while (fgets (textLine, TID_MAXTEXTSTRING, filePtr))
    {
        stringLength = strlen (textLine);

        // check to see if there's a space after the last data entry on the line. if so, our space counter loop below should stop prior to that final space. this allows us to read text files written both with and without spaces after the final entry of a line
        // stringLength-2 would be the position for this space, because the final character is a carriage return (10) at position stringLength - 1
        if (textLine[stringLength - 2] == 32)
            stringLength = stringLength - 2; // lop off the final space and CR

        numSpaces = 0;

        for (j = 0; j < stringLength; j++)
        {
            if (textLine[j] == 32)
                numSpaces++;
        }

        // there's a space after each entry in a file written by write_text(), except for the final entry. So (numSpaces + 1)==length
        x->x_instances[i].length = numSpaces + 1;

        if (x->x_instances[i].length > maxLength)
            maxLength = x->x_instances[i].length;

        if (x->x_instances[i].length < minLength)
            minLength = x->x_instances[i].length;

        // get the appropriate number of bytes for the data
        x->x_instances[i].data = (t_float *)t_getbytes (x->x_instances[i].length * sizeof (t_float));

        i++;
    }

    x->x_minFeatureLength = minLength;
    x->x_maxFeatureLength = maxLength;
    x->x_neighborhood = x->x_numInstances;
    x->x_numClusters = x->x_numInstances;

    // update x_attributeData based on new x_maxFeatureLength. postFlag argument TRUE here
    timbreID_attributeDataResize (x, 0, x->x_maxFeatureLength, 1);

    for (i = 0; i < x->x_numInstances; i++)
    {
        x->x_clusters[i].numMembers = 2;
        x->x_clusters[i].members = (t_instanceIdx *)t_getbytes (x->x_clusters[i].numMembers * sizeof (t_instanceIdx));

        x->x_clusters[i].members[0] = i; // first member of the cluster is the instance index
        x->x_clusters[i].members[1] = UINT_MAX;

        x->x_instances[i].clusterMembership = i; // init instance's cluster membership to index
    }

    // now that we have the number of instances and the length of each feature, close and re-open file to extract the actual data
    fclose (filePtr);
    filePtr = fopen (fileNameBuf, "r");

    for (i = 0; i < x->x_numInstances; i++)
    {
        featurePtr = x->x_instances[i].data;

        for (j = 0; j < x->x_instances[i].length; j++, featurePtr++)
            fscanf (filePtr, "%f", featurePtr);
    };

    post ("%s: read %i instances from %s.\n", x->x_objSymbol->s_name, x->x_numInstances, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_ARFF (t_timbreID* x, t_symbol* s, int argc, t_atom* argv)
{
    FILE* filePtr;
    t_instanceIdx i, j, attRangeLow, attRangeHi, argCount;
    t_float* featurePtr;
    t_symbol* filenameSymbol;
    t_symbol* relationSymbol;
    t_symbol* attSymbol;
    const char* fileName;
    const char* relation;
    const char* attName;
    char fileNameBuf[MAXPDSTRING];

    attRangeLow = 0;
    attRangeHi = 0;
    attSymbol = 0;
    attName = 0;

    argCount = (argc < 0) ? 0 : argc;

    filenameSymbol = atom_getsymbol (argv);
    fileName = filenameSymbol->s_name;

    canvas_makefilename (x->x_canvas, fileName, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "w");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    if (argCount > 1)
    {
        relationSymbol = atom_getsymbol (argv + 1);
        relation = relationSymbol->s_name;
    }
    else
        relation = "relation";

    fprintf (filePtr, "@RELATION %s\n\n\n", relation);

    if (argCount > 2)
    {
        for (i = 2; i < argCount; i++)
        {

            switch ((i - 2) % 3)
            {
                case 0:
                    attRangeLow = atom_getfloat (argv + i);
                    break;
                case 1:
                    attRangeHi = atom_getfloat (argv + i);
                    break;
                case 2:
                    attSymbol = atom_getsymbol (argv + i);
                    attName = attSymbol->s_name;

                    for (j = 0; j <= attRangeHi - attRangeLow; j++)
                        fprintf (filePtr, "@ATTRIBUTE %s-%i NUMERIC\n", attName, j);
                    break;
                default:
                    break;
            }
        }

        // BUG: this was causing a crash because x_maxFeatureLength and attRangeHi are unsigned integers, so we won't get a negative number as expected when there are indeed enough arguments. This came up with version 0.7 because of typedefs - no longer using int. Quick fix is to typecast back to int during the arithmetic

        // in case the argument list was incomplete
        if ((x->x_maxFeatureLength - 1) > attRangeHi)
        {
            for (i = attRangeHi + 1, j = 0; i < x->x_maxFeatureLength; i++, j++)
                fprintf (filePtr, "@ATTRIBUTE undefined-attribute-%i NUMERIC\n", j);
        }
    }
    else
    {
        for (i = 0; i < x->x_maxFeatureLength; i++)
            fprintf (filePtr, "@ATTRIBUTE undefined-attribute-%i NUMERIC\n", i);
    }

    fprintf (filePtr, "\n\n");
    fprintf (filePtr, "@DATA\n\n");

    for (i = 0; i < x->x_numInstances; i++)
    {
        featurePtr = x->x_instances[i].data;

        for (j = 0; j < x->x_instances[i].length - 1; j++)
            fprintf (filePtr, "%0.20f, ", *featurePtr++);

        // last attribute shouldn't be followed by a comma and space
        fprintf (filePtr, "%0.20f", *featurePtr++);

        fprintf (filePtr, "\n");
       };

    post ("%s: wrote %i non-normalized instances to %s.", x->x_objSymbol->s_name, x->x_numInstances, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_MATLAB (t_timbreID* x, t_symbol* file_symbol, t_symbol* var_symbol)
{
    FILE* filePtr;
    t_instanceIdx i, featuresWritten;
    t_float* featurePtr;
    char fileNameBuf[MAXPDSTRING];

    canvas_makefilename (x->x_canvas, file_symbol->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "w");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    if (x->x_maxFeatureLength != x->x_minFeatureLength)
    {
        pd_error (x, "%s: database instances must have uniform length for MATLAB matrix formatting. failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        fclose (filePtr);
        return;
    }

    fprintf (filePtr, "%% name: %s\n", var_symbol->s_name);
    fprintf (filePtr, "%% type: matrix\n");
    fprintf (filePtr, "%% rows: %i\n", x->x_numInstances);
    fprintf (filePtr, "%% columns: %i\n\n", x->x_maxFeatureLength);

    for (i = 0; i < x->x_numInstances; i++)
    {
        featurePtr = x->x_instances[i].data;

        featuresWritten = 0; // to keep track of each instances no. of features written.

        while (1)
        {
            if (featuresWritten++ == x->x_instances[i].length - 1)
            {
                fprintf (filePtr, "%0.20f", *featurePtr++);
                break;
            }
            else
                fprintf (filePtr, "%0.20f, ", *featurePtr++);
        };

        fprintf (filePtr, "\n");
       };

    post ("%s: wrote %i non-normalized instances to %s.", x->x_objSymbol->s_name, x->x_numInstances, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_OCTAVE (t_timbreID* x, t_symbol* file_symbol, t_symbol* var_symbol)
{
    FILE* filePtr;
    t_instanceIdx i, featuresWritten;
    t_float* featurePtr;
    char fileNameBuf[MAXPDSTRING];

    canvas_makefilename (x->x_canvas, file_symbol->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "w");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    if (x->x_maxFeatureLength != x->x_minFeatureLength)
    {
        pd_error (x, "%s: database instances must have uniform length for MATLAB matrix formatting. failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        fclose (filePtr);
        return;
    }

    fprintf (filePtr, "# Created by timbreID version %s\n", TID_VERSION);
    fprintf (filePtr, "# name: %s\n", var_symbol->s_name);
    fprintf (filePtr, "# type: matrix\n");
    fprintf (filePtr, "# rows: %i\n", x->x_numInstances);
    fprintf (filePtr, "# columns: %i\n", x->x_maxFeatureLength);

    for (i = 0; i < x->x_numInstances; i++)
    {
        featurePtr = x->x_instances[i].data;

        featuresWritten = 0; // to keep track of each instances no. of features written.

        while (1)
        {
            if (featuresWritten++ == x->x_instances[i].length - 1)
            {
                fprintf (filePtr, " %0.20f\n", *featurePtr++);
                break;
            }
            else
                fprintf (filePtr, " %0.20f", *featurePtr++);
        };
    };

    post ("%s: wrote %i non-normalized instances to %s.", x->x_objSymbol->s_name, x->x_numInstances, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_FANN (t_timbreID* x, t_symbol* s, t_floatarg normRange)
{
    FILE* filePtr;
    t_instanceIdx i, j;
    char fileNameBuf[MAXPDSTRING];

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "w");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    // write the header containing the number of instances, input size, output size
    fprintf (filePtr, "%i ", x->x_numInstances);
    // note that this is minFeatureLength assuming that all features will be of the same length so x_minFeatureLength==x_maxFeatureLength
    fprintf (filePtr, "%i ", x->x_minFeatureLength);
    fprintf (filePtr, "%i ", x->x_numClusters);
    fprintf (filePtr, "\n");

    normRange = (normRange < 0) ? 0 : normRange;
    normRange = (normRange > 1) ? 1 : normRange;

    for (i = 0; i < x->x_numInstances; i++)
    {
        for (j = 0; j < x->x_minFeatureLength; j++)
        {
            if (x->x_normalize)
            {
                if (normRange)
                    fprintf (filePtr, "%f ", ((x->x_instances[i].data[j] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar * 2.0) - 1.0);
                else
                    fprintf (filePtr, "%f ", (x->x_instances[i].data[j] - x->x_attributeData[j].normData.minVal) * x->x_attributeData[j].normData.normScalar);
            }
            else
                fprintf (filePtr, "%f ", (x->x_instances[i].data[j]));
        }

        // carriage return between input data line and label line
        fprintf (filePtr, "\n");

        // TODO: fill with all zeros, except a 1.0 in the slot of the cluster idx
        for (j = 0; j < x->x_numClusters; j++)
        {
            if (j == x->x_instances[i].clusterMembership)
                fprintf (filePtr, "%i ", 1);
            else
                fprintf (filePtr, "%i ", 0);
        }

        // carriage return between label line and next line of data
        fprintf (filePtr, "\n");
       };

    post ("%s: wrote %i training instances and labels to %s.", x->x_objSymbol->s_name, x->x_numInstances, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_writeClusters (t_timbreID* x, t_symbol* s)
{
    FILE* filePtr;
    t_instanceIdx i;
    char fileNameBuf[MAXPDSTRING];

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "wb");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    // write a header indicating the number of clusters (x->x_numClusters)
    fwrite (&x->x_numClusters, sizeof (t_instanceIdx), 1, filePtr);

    // write the cluster memberships next
    for (i = 0; i < x->x_numInstances; i++)
        fwrite (&x->x_instances[i].clusterMembership, sizeof (t_instanceIdx), 1, filePtr);

    // write the number of members for each cluster
    for (i = 0; i < x->x_numInstances; i++)
        fwrite (&x->x_clusters[i].numMembers, sizeof (t_instanceIdx), 1, filePtr);

    // finally, write the members data
    for (i = 0; i < x->x_numInstances; i++)
        fwrite (x->x_clusters[i].members, sizeof (t_instanceIdx), x->x_clusters[i].numMembers, filePtr);

    post ("%s: wrote %i clusters to %s.", x->x_objSymbol->s_name, x->x_numClusters, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_readClusters (t_timbreID* x, t_symbol* s)
{
    FILE* filePtr;
    t_instanceIdx i, numClusters;
    char fileNameBuf[MAXPDSTRING];
    int cTest;
    t_bool txtFlag;

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "rb");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to open %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    // test to see if it's a .clu or .txt file. if any fgetc() results are > 127, it's not a text file.
    txtFlag = true;

    while ((cTest = fgetc (filePtr)) != EOF)
    {
        if (cTest > 127)
            txtFlag = false;
    }

    if (txtFlag)
    {
        pd_error (x, "%s: data in file %s is plain text (.txt), not binary (.clu). use read_clusters_text instead.", x->x_objSymbol->s_name, fileNameBuf);
        fclose (filePtr);
        // TODO: shouldn't have an else branch after a return statement
        return;
    }
    else
    {
        // if we're safe, close and re-open so we can get the binary data.
        fclose (filePtr);
        filePtr = fopen (fileNameBuf, "rb");
    }

    // read header indicating number of clusters
    fread (&numClusters, sizeof (t_instanceIdx), 1, filePtr);

    if (numClusters > x->x_numInstances)
    {
        pd_error (x, "%s: %s contains more clusters than current number of database instances. read failed.", x->x_objSymbol->s_name, fileNameBuf);
        fclose (filePtr);
        return;
    }

    x->x_numClusters = numClusters;

    // read the cluster memberships next
    for (i = 0; i < x->x_numInstances; i++)
        fread (&x->x_instances[i].clusterMembership, sizeof (t_instanceIdx), 1, filePtr);

    // free members memory, read the number of members and get new members memory
    for (i = 0; i < x->x_numInstances; i++)
    {
        t_freebytes (x->x_clusters[i].members, x->x_clusters[i].numMembers * sizeof (t_instanceIdx));
        fread (&x->x_clusters[i].numMembers, sizeof (t_instanceIdx), 1, filePtr);
        x->x_clusters[i].members = (t_instanceIdx *)t_getbytes (x->x_clusters[i].numMembers * sizeof (t_instanceIdx));
    }

    // read the actual members data
    for (i = 0; i < x->x_numInstances; i++)
        fread (x->x_clusters[i].members, sizeof (t_instanceIdx), x->x_clusters[i].numMembers, filePtr);

    post ("%s: read %i clusters from %s.", x->x_objSymbol->s_name, x->x_numClusters, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_writeClustersText (t_timbreID* x, t_symbol* s)
{
    FILE* filePtr;
    t_instanceIdx i, j;
    t_instanceIdx* clusterPtr;
    char fileNameBuf[MAXPDSTRING];

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "w");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to create %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    for (i = 0; i < x->x_numClusters; i++)
    {
        clusterPtr = x->x_clusters[i].members;

        for (j = 0; j < x->x_clusters[i].numMembers-2; j++)
            fprintf (filePtr, "%i ", *clusterPtr++);

        // no space for the final instance given on the line for cluster i
        fprintf (filePtr, "%i", *clusterPtr++);

        // newline to end the list of instances for cluster i
        fprintf (filePtr, "\n");
    };

    post ("%s: wrote %i clusters to %s.", x->x_objSymbol->s_name, x->x_numClusters, fileNameBuf);

    fclose (filePtr);
}


static void timbreID_readClustersText (t_timbreID* x, t_symbol* s)
{
    FILE* filePtr;
    t_instanceIdx i, j, stringLength, numSpaces, numClusters;
    t_instanceIdx* featurePtr;
    char fileNameBuf[MAXPDSTRING], textLine[TID_MAXTEXTSTRING];

    canvas_makefilename (x->x_canvas, s->s_name, fileNameBuf, MAXPDSTRING);

    filePtr = fopen (fileNameBuf, "r");

    if ( !filePtr)
    {
        pd_error (x, "%s: failed to open %s", x->x_objSymbol->s_name, fileNameBuf);
        return;
    }

    // count lines in text file to determine numClusters
    numClusters = 0;

    while (fgets (textLine, TID_MAXTEXTSTRING, filePtr))
        numClusters++;

    if (numClusters > x->x_numInstances)
    {
        pd_error (x, "%s: %s contains more clusters than current number of database instances. read failed.", x->x_objSymbol->s_name, fileNameBuf);
        fclose (filePtr);
        return;
    }

    x->x_numClusters = numClusters;

    // free existing x->x_clusters[i].members memory
    for (i = 0; i < x->x_numInstances; i++)
        t_freebytes (x->x_clusters[i].members, x->x_clusters[i].numMembers * sizeof (t_instanceIdx));


    fclose (filePtr);
    filePtr = fopen (fileNameBuf, "r");

    i = 0;

    while (fgets (textLine, TID_MAXTEXTSTRING, filePtr))
    {
        stringLength = strlen (textLine);

        // check to see if there's a space after the last data entry on the line. if so, our space counter loop below should stop prior to that final space. this allows us to read text files written both with and without spaces after the final entry of a line
        // stringLength-2 would be the position for this space, because the final character is a carriage return (10) at position stringLength - 1
        if (textLine[stringLength - 2] == 32)
            stringLength = stringLength - 2; // lop off the final space and CR

        numSpaces = 0;

        for (j = 0; j < stringLength; j++)
        {
            if (textLine[j] == 32)
                numSpaces++;
        }

        // there's a space after each entry in a file written by write_clusters_text(), except for the final instance on the line. So numMembers should be numSpaces + 1. However, we must add one more slot for the terminating UINT_MAX of each cluster list. So in the end, numMembers should be numSpaces+2
        x->x_clusters[i].numMembers = numSpaces + 2;

        // get the appropriate number of bytes for the data
       x->x_clusters[i].members = (t_instanceIdx *)t_getbytes (x->x_clusters[i].numMembers * sizeof (t_instanceIdx));

        i++;
    }

    fclose (filePtr);
    filePtr = fopen (fileNameBuf, "r");

    for (i = 0; i < x->x_numClusters; i++)
    {
        featurePtr = x->x_clusters[i].members;

        for (j = 0; j < x->x_clusters[i].numMembers - 1; j++, featurePtr++)
            fscanf (filePtr, "%i", featurePtr);

        // don't forget to terminate with UINT_MAX
        x->x_clusters[i].members[j] = UINT_MAX;
    };

    // fill out any remaining with the default setup
    for (; i < x->x_numInstances; i++)
    {
        x->x_clusters[i].numMembers = 2;
        x->x_clusters[i].members = (t_instanceIdx *)t_getbytes (x->x_clusters[i].numMembers * sizeof (t_instanceIdx));
        x->x_clusters[i].members[0] = i;
        x->x_clusters[i].members[1] = UINT_MAX;

        // set cluster membership to instance index
        x->x_instances[i].clusterMembership = i;
    }

    for (i = 0; i < x->x_numClusters; i++)
    {
        j = 0;

        while (x->x_clusters[i].members[j] != UINT_MAX)
        {
            t_instanceIdx idx;
            idx = x->x_clusters[i].members[j];
            x->x_instances[idx].clusterMembership = i;
            j++;
        };
    };

    post ("%s: read %i clusters from %s", x->x_objSymbol->s_name, x->x_numClusters, fileNameBuf);

    fclose (filePtr);
}


static void* timbreID_new (void)
{
    t_timbreID* x = (t_timbreID *)pd_new (timbreID_class);

    // define outlets
    x->x_id = outlet_new (&x->x_obj, &s_float);
    x->x_nearestDist = outlet_new (&x->x_obj, &s_float);
    x->x_confidence = outlet_new (&x->x_obj, &s_float);
    x->x_listOut = outlet_new (&x->x_obj, gensym ("list"));
    inlet_new (&x->x_obj, &x->x_obj.ob_pd, gensym ("list"), gensym ("id"));
    inlet_new (&x->x_obj, &x->x_obj.ob_pd, gensym ("list"), gensym ("concat_id"));

    x->x_objSymbol = gensym ("timbreID");

    x->x_instances = (t_instance *)t_getbytes (0);
    x->x_classRefs = (t_instance *)t_getbytes (0);
    x->x_clusters = (t_cluster *)t_getbytes (0);
    x->x_attributeData = (t_attributeData *)t_getbytes (0);

    x->x_maxFeatureLength = 0;
    x->x_minFeatureLength = INT_MAX;
    x->x_numClusters=0;
    x->x_numInstances = 0;
    x->x_numClassRefs = 0;
    x->x_distMetric = euclidean;
    x->x_k = 1;
    x->x_outputKnnMatches = false;
    x->x_normalize = false;
    x->x_relativeOrdering = true;
    x->x_stutterProtect = false;

    x->x_concatWrap = true;
    x->x_prevMatch = UINT_MAX;
    x->x_maxMatches = 3;
    x->x_reorientFlag = false;
    x->x_neighborhood = 0;
    x->x_searchCenter = 0;
    x->x_jumpProb = 0.0;

    x->x_canvas = canvas_getcurrent();

    return (x);
}


static void timbreID_free (t_timbreID* x)
{
    t_instanceIdx i;

    // free the database memory per instance
    for (i = 0; i < x->x_numInstances; i++)
        t_freebytes (x->x_instances[i].data, x->x_instances[i].length * sizeof (t_float));

    t_freebytes (x->x_instances, x->x_numInstances * sizeof (t_instance));

    // free the class reference memory per instance
    for (i = 0; i < x->x_numClassRefs; i++)
        t_freebytes (x->x_classRefs[i].data, x->x_classRefs[i].length * sizeof (t_float));

    t_freebytes (x->x_classRefs, x->x_numClassRefs * sizeof (t_instance));

    // free the cluster memory
    for (i = 0; i < x->x_numInstances; i++)
        t_freebytes (x->x_clusters[i].members, x->x_clusters[i].numMembers * sizeof (t_instanceIdx));

    t_freebytes (x->x_clusters, x->x_numInstances * sizeof (t_cluster));

    // free the attribute data memory
    t_freebytes (x->x_attributeData, x->x_maxFeatureLength * sizeof (t_attributeData));
}


void timbreID_setup (void)
{
    timbreID_class =
    class_new (
        gensym ("timbreID"),
        (t_newmethod)timbreID_new,
        (t_method)timbreID_free,
        sizeof (t_timbreID),
        CLASS_DEFAULT,
        0
    );

    class_addcreator (
        (t_newmethod)timbreID_new,
        gensym ("timbreIDLib/timbreID"),
        0
    );

    class_addlist (
        timbreID_class,
        (t_method)timbreID_train
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_id,
        gensym ("id"),
        A_GIMME,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_worstMatch,
        gensym ("worst_match"),
        A_GIMME,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_concatId,
        gensym ("concat_id"),
        A_GIMME,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_concatSearchWrap,
        gensym ("concat_search_wrap"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_concatNeighborhood,
        gensym ("neighborhood"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_concatJumpProb,
        gensym ("jump_prob"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_concatReorient,
        gensym ("reorient"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_concatSearchCenter,
        gensym ("search_center"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_concatMaxMatches,
        gensym ("max_matches"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_concatStutterProtect,
        gensym ("stutter_protect"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_knn,
        gensym ("knn"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_outputKnnMatches,
        gensym ("output_knn_matches"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_normalize,
        gensym ("normalize"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_print,
        gensym ("print"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_numInstances,
        gensym ("num_instances"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_manualCluster,
        gensym ("manual_cluster"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_computeCluster,
        gensym ("cluster"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_uncluster,
        gensym ("uncluster"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_clustersList,
        gensym ("clusters_list"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_clusterList,
        gensym ("cluster_list"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_clusterMembership,
        gensym ("cluster_membership"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_computeOrder,
        gensym ("order"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_relativeOrdering,
        gensym ("relative_ordering"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_computeVariance,
        gensym ("order_attributes_by_variance"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_distMetric,
        gensym ("dist_metric"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_weights,
        gensym ("weights"),
        A_GIMME,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_attributeNames,
        gensym ("attribute_names"),
        A_GIMME,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_attributeInfo,
        gensym ("attribute_info"),
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_attributeOrder,
        gensym ("attribute_order"),
        A_GIMME,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_attributeRange,
        gensym ("attribute_range"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_reorderAttributes,
        gensym ("reorder_attributes"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_featureList,
        gensym ("feature_list"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

// adding this additional method name to retrieve an instance. makes it more consistent with new attribute_list name. adding a warning that feature_list is deprecated
    class_addmethod (
        timbreID_class,
        (t_method)timbreID_instanceList,
        gensym ("instance_list"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_attributeList,
        gensym ("attribute_list"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_similarityMatrix,
        gensym ("similarity_matrix"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_maxValues,
        gensym ("max_values"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_minValues,
        gensym ("min_values"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_maxk,
        gensym ("max_k_values"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_mink,
        gensym ("min_k_values"),
        A_DEFFLOAT,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_classReference,
        gensym ("class_reference"),
        A_GIMME,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_maxFeatureLength,
        gensym ("max_feature_length"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_minFeatureLength,
        gensym ("min_feature_length"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_clear,
        gensym ("clear"),
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_write,
        gensym ("write"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_read,
        gensym ("read"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_writeText,
        gensym ("write_text"),
        A_SYMBOL,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_readText,
        gensym ("read_text"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_ARFF,
        gensym ("ARFF"),
        A_GIMME,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_MATLAB,
        gensym ("MATLAB"),
        A_SYMBOL,
        A_SYMBOL,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_OCTAVE,
        gensym ("OCTAVE"),
        A_SYMBOL,
        A_SYMBOL,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_FANN,
        gensym ("FANN"),
        A_SYMBOL,
        A_DEFFLOAT,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_writeClusters,
        gensym ("write_clusters"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_readClusters,
        gensym ("read_clusters"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_writeClustersText,
        gensym ("write_clusters_text"),
        A_SYMBOL,
        0
    );

    class_addmethod (
        timbreID_class,
        (t_method)timbreID_readClustersText,
        gensym ("read_clusters_text"),
        A_SYMBOL,
        0
    );
}
