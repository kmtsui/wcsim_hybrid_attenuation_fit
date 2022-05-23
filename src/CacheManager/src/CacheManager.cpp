#include "CacheManager.h"

#include <vector>
#include <set>


Cache::Manager* Cache::Manager::fSingleton = nullptr;
std::map<const AnaFitParameters*, int> Cache::Manager::ParameterMap;

Cache::Manager::Manager(int events, int parameters,
                        int norms,
                        int attens,
                        int attenzs,
                        int polys,
                        int sphis,
                        int splines, int splinePoints,
                        int histBins) {
    std::cout << "Creating cache manager" << std::endl;

    fTotalBytes = 0;
    try {
        fParameterCache.reset(new Cache::Parameters(parameters));
        fTotalBytes += fParameterCache->GetResidentMemory();

        fWeightsCache.reset(
            new Cache::Weights(fParameterCache->GetParameters(),
                               fParameterCache->GetLowerClamps(),
                               fParameterCache->GetUpperClamps(),
                               events));
        fTotalBytes += fWeightsCache->GetResidentMemory();

        fNormalizations.reset(new Cache::Weight::Normalization(
                                  fWeightsCache->GetWeights(),
                                  fParameterCache->GetParameters(),
                                  norms));
        fWeightsCache->AddWeightCalculator(fNormalizations.get());
        fTotalBytes += fNormalizations->GetResidentMemory();

        fAttenuations.reset(new Cache::Weight::Attenuation(
                                fWeightsCache->GetWeights(),
                                fParameterCache->GetParameters(),
                                attens));
        fWeightsCache->AddWeightCalculator(fAttenuations.get());
        fTotalBytes += fAttenuations->GetResidentMemory();

        fAttenuationzs.reset(new Cache::Weight::AttenuationZ(
                                fWeightsCache->GetWeights(),
                                fParameterCache->GetParameters(),
                                attenzs));
        fWeightsCache->AddWeightCalculator(fAttenuationzs.get());
        fTotalBytes += fAttenuationzs->GetResidentMemory();

        fPolynomialCosth.reset(new Cache::Weight::PolynomialCosth(
                                   fWeightsCache->GetWeights(),
                                   fParameterCache->GetParameters(),
                                   polys));
        fWeightsCache->AddWeightCalculator(fPolynomialCosth.get());
        fTotalBytes += fPolynomialCosth->GetResidentMemory();

        fSourcePhiVars.reset(new Cache::Weight::SourcePhiVar(
                                   fWeightsCache->GetWeights(),
                                   fParameterCache->GetParameters(),
                                   sphis));
        fWeightsCache->AddWeightCalculator(fSourcePhiVars.get());
        fTotalBytes += fSourcePhiVars->GetResidentMemory();

        fSplines.reset(new Cache::Weight::Spline(
                                  fWeightsCache->GetWeights(),
                                  fParameterCache->GetParameters(),
                                  fParameterCache->GetLowerClamps(),
                                  fParameterCache->GetUpperClamps(),
                                  splines, splinePoints));
        fWeightsCache->AddWeightCalculator(fSplines.get());
        fTotalBytes += fSplines->GetResidentMemory();

        fHistogramsCache.reset(new Cache::IndexedSums(
                                  fWeightsCache->GetWeights(),
                                  histBins));
        fTotalBytes += fHistogramsCache->GetResidentMemory();

    }
    catch (...) {
        std::cout << "Failed to allocate memory, so stopping" << std::endl;
        throw std::runtime_error("Not enough memory available");
    }

    std::cout << "Approximate cache manager size for"
            << " " << events << " events:"
            << " " << GetResidentMemory()/1E+9 << " GB "
            << " (" << GetResidentMemory()/events << " bytes per event)"
            << std::endl;
}

bool Cache::Manager::Build(std::vector<AnaSample*> samples, std::vector<AnaFitParameters*> fitparas) {
    std::cout << "Build the cache for Cache::Manager" << std::endl;

    int events = 0;
    int norms = 0;
    int attens = 0;
    int attenzs = 0;
    int polys = 0;
    int sphis = 0;
    int splines = 0;
    int splinePoints = 0;
    Cache::Manager::ParameterMap.clear();

    std::set<const AnaFitParameters*> usedParameters;
    for(auto& sample : samples){
        std::cout << "Sample " << sample->GetName()
                << " with " << sample->GetNPMTs()
                << " events" << std::endl;
        events += sample->GetNPMTs();
        for (int i=0; i<sample->GetNPMTs(); i++){
            // The reduce index to save the result for this event.
            AnaEvent* event = sample->GetPMT(i);
            if (sample->UseTemplate()) events += event->GetTimetofNom().size();
            for(auto& fitpara : fitparas)
            {
                if ( fitpara->GetParameterFunctionType() == kIdentity )
                {
                    int bin = fitpara->GetParBin(event->GetSampleType(), i);
                    if(bin == PASSEVENT || bin == BADBIN) continue;
                    norms++;
                }
                else if ( fitpara->GetParameterFunctionType() == kAttenuation )
                {
                    int bin = fitpara->GetParBin(event->GetSampleType(), i);
                    if(bin == PASSEVENT || bin == BADBIN) continue;
                    attens++;
                }
                else if ( fitpara->GetParameterFunctionType() == kAttenuationZ )
                {
                    int bin = fitpara->GetParBin(event->GetSampleType(), i);
                    if(bin == PASSEVENT || bin == BADBIN) continue;
                    attenzs++;
                }
                else if ( fitpara->GetParameterFunctionType() == kPolynomialCosth )
                {
                    int bin = fitpara->GetParBin(event->GetSampleType(), i);
                    if(bin == PASSEVENT || bin == BADBIN) continue;
                    polys++;
                }
                else if ( fitpara->GetParameterFunctionType() == kSourcePhiVar )
                {
                    int bin = fitpara->GetParBin(event->GetSampleType(), i);
                    if(bin == PASSEVENT || bin == BADBIN) continue;
                    sphis++;
                }

                if ( fitpara->UseSpline() && sample->UseTemplate())
                {
                    std::vector<TGraph*> graphs = fitpara->GetSplineGraph(event->GetSampleType(), i);
                    for(auto& graph : graphs)
                    {
                        splines++;
                        splinePoints += graph->GetN();
                    }
                }
            }
        }
    }

    // Count the total number of histogram cells.
    int histCells = 0;
    for(auto& sample : samples){
        int cells = sample->GetNBins();
        std::cout << "Add histogram for " << sample->GetName()
                << " with " << cells
                << " cells" << std::endl;
        histCells += cells;

        if (sample->UseTemplate())
        {
            int ntcells = sample->GetNTemplateBins();
            std::cout << "Add template histogram for " << sample->GetName()
                    << " with " << ntcells
                    << " cells" << std::endl;
            histCells += ntcells;
        }
    }

    //int parameters = fitparas.size();
    int parameters = 0;
    for(auto& fitpara : fitparas)
    {
        if ( fitpara->GetParameterFunctionType() == kPolynomialCosth )
        {
            std::vector<int> pol_orders = ((PolynomialCosth*)fitpara->GetParameterFunction())->pol_orders;
            for (int i=0;i<pol_orders.size();i++)   
                parameters += pol_orders[i]+1;
        }
        else 
            parameters += fitpara->GetNpar();
    }
    std::cout << "Cache for " << events << " events --"
            << " using " << parameters << " parameters"
            << std::endl;
    std::cout << "    Histogram bins: " << histCells
            << " (" << 1.0*events/histCells << " events per bin)"
            << std::endl;

    // Try to allocate the GPU
    if (!Cache::Manager::Get()) {
        fSingleton = new Manager(events,parameters,
                                 norms,
                                 attens,
                                 attenzs,
                                 polys,
                                 sphis,
                                 splines, splinePoints,
                                 histCells);
    }

    // In case the cache isn't allocated (usually because it's turned off on
    // the command line).
    if (!Cache::Manager::Get()) {
        std::cout << "Cache will not be used"
                << std::endl;
        return false;
    }

    // Add the dials to the cache.
    int usedResults = 0; // Number of cached results that have been used up.
    int parCount = 0;
    for(auto& fitpara : fitparas)
    {
        Cache::Manager::ParameterMap[fitpara] = parCount;

        if ( fitpara->GetParameterFunctionType() == kPolynomialCosth )
        {
            std::vector<int> pol_orders = ((PolynomialCosth*)fitpara->GetParameterFunction())->pol_orders;
            for (int i=0;i<pol_orders.size();i++)   
                parCount += pol_orders[i]+1;
        }
        else
            parCount += fitpara->GetNpar();
    }
    for(auto& sample : samples) {
        std::cout << "Fill cache for " << sample->GetName()
                << " with " << sample->GetNPMTs()
                << " events" << std::endl;
        for (int i=0; i<sample->GetNPMTs(); i++)
        {
            // The reduce index to save the result for this event.
            AnaEvent* event = sample->GetPMT(i);
            int resultIndex = usedResults++;
            event->setCacheManagerIndex(resultIndex);
            event->setCacheManagerValuePointer(Cache::Manager::Get()
                                              ->GetWeightsCache()
                                              .GetResultPointer(resultIndex));
            event->setCacheManagerValidPointer(Cache::Manager::Get()
                                              ->GetWeightsCache()
                                              .GetResultValidPointer());
            event->setCacheManagerUpdatePointer(
                [](){Cache::Manager::Get()->GetWeightsCache().GetResult(0);});
            Cache::Manager::Get()
                ->GetWeightsCache()
                .SetInitialValue(resultIndex,event->GetEvWghtMC());
            if (sample->UseTemplate()) 
            {
                std::vector<double> timetofNom = event->GetTimetofNom();
                usedResults += timetofNom.size();
                for (int j=0;j<timetofNom.size();j++)
                {
                    Cache::Manager::Get()
                        ->GetWeightsCache()
                        .SetInitialValue(resultIndex+j+1,timetofNom[j]);
                }
            }
            for(auto& fitpara : fitparas)
            {
                int parIndex = Cache::Manager::ParameterMap[fitpara];

		        int bin = fitpara->GetParBin(event->GetSampleType(), i);
                if(bin == PASSEVENT || bin == BADBIN) continue;

                if ( fitpara->GetParameterFunctionType() == kIdentity )
                {
                    Cache::Manager::Get()
                        ->fNormalizations
                        ->ReserveNorm(resultIndex,parIndex+bin);
                }
                else if ( fitpara->GetParameterFunctionType() == kAttenuation )
                {
                    Cache::Manager::Get()
                        ->fAttenuations
                        ->ReserveAtten(resultIndex,parIndex+bin,event->GetR(),event->GetOmega()*event->GetEff());
                }
                else if ( fitpara->GetParameterFunctionType() == kAttenuationZ )
                {
                    Cache::Manager::Get()
                        ->fAttenuationzs
                        ->ReserveAttenz(resultIndex,parIndex+bin,event->GetR(),event->GetOmega()*event->GetEff(),event->GetZ0(),event->GetDz());
                }
                else if ( fitpara->GetParameterFunctionType() == kPolynomialCosth )
                {
                    std::vector<int>    pol_orders = ((PolynomialCosth*)fitpara->GetParameterFunction())->pol_orders;
                    std::vector<double> pol_range = ((PolynomialCosth*)fitpara->GetParameterFunction())->pol_range;
                    double costh = event->GetCosth();
                    int polyOrd = pol_orders[0];
                    double polyArg = costh-pol_range[0];
                    if (polyArg>0)
                    { 
                        for (int j=0;j<pol_orders.size();j++)
                        {
                            if ((costh>=pol_range[j] && costh<pol_range[j+1]) || j==pol_orders.size()-1)
                            {
                                polyArg = costh-pol_range[j];
                                polyOrd = pol_orders[j];
                                break;
                            }
                            parIndex += pol_orders[j] + 1;
                        }
                    }
                    Cache::Manager::Get()
                        ->fPolynomialCosth
                        ->ReservePoly(resultIndex,parIndex,polyOrd,polyArg);
                }
                else if ( fitpara->GetParameterFunctionType() == kSourcePhiVar )
                {
                    Cache::Manager::Get()
                        ->fSourcePhiVars
                        ->ReserveSPhi(resultIndex,parIndex+bin,cos(event->GetPhis()));
                }

                if ( fitpara->UseSpline() && sample->UseTemplate())
                {
                    std::vector<TGraph*> graphs = fitpara->GetSplineGraph(event->GetSampleType(), i);
                    for(int j=0;j<graphs.size();j++)
                    {
                        Cache::Manager::Get()
                        ->fSplines
                        ->AddSpline(resultIndex+j+1,parIndex+bin,graphs[j]);
                    }
                }
            }
        }
    }

    // Error checking!
    if (usedResults != Cache::Manager::Get()
        ->GetWeightsCache().GetResultCount()) {
        std::cout << "Cache Manager -- used Results:     "
                 << usedResults << std::endl;
        std::cout << "Cache Manager -- expected Results: "
                 << Cache::Manager::Get()->GetWeightsCache().GetResultCount()
                 << std::endl;
        throw std::runtime_error("Probable problem putting dials in cache");
    }

    // Add this histogram cells to the cache.
    int nextHist = 0;
    for(auto& sample : samples) {
        std::cout << "Fill cache for " << sample->GetName()
                << " with " << sample->GetNPMTs()
                << " events" << std::endl;
        int thisHist = nextHist;
        sample->setCacheManagerIndex(thisHist);
        sample->setCacheManagerValuePointer(
            Cache::Manager::Get()->GetHistogramsCache()
            .GetSumsPointer());
        sample->setCacheManagerValidPointer(
            Cache::Manager::Get()->GetHistogramsCache()
            .GetSumsValidPointer());
        sample->setCacheManagerUpdatePointer(
            [](){Cache::Manager::Get()->GetHistogramsCache().GetSum(0);});
        int cells = sample->GetNBins();
        nextHist += cells;
        for (int i=0; i<sample->GetNPMTs(); i++){
            AnaEvent* event = sample->GetPMT(i);
            int eventIndex = event->getCacheManagerIndex();
            int cellIndex = event->GetSampleBin();
            if (cellIndex < 0 || cells <= cellIndex) {
                throw std::runtime_error("Histogram bin out of range");
            }
            int theEntry = thisHist + cellIndex;
            Cache::Manager::Get()->GetHistogramsCache()
                .SetEventIndex(eventIndex,theEntry);
        }

        if (sample->UseTemplate())
        {
            thisHist = nextHist;
            int ntcells = sample->GetNTemplateBins();
            nextHist += ntcells;
            for (int i=0; i<sample->GetNPMTs(); i++)
            {
                AnaEvent* event = sample->GetPMT(i);
                int eventIndex = event->getCacheManagerIndex();
                for (int j=0;j<sample->GetTemplate()->GetNbinsY();j++)
                {
                    int cellIndex =  sample->CombineTemplate() ? j : cells*j + event->GetSampleBin();
                    if (cellIndex < 0 || ntcells <= cellIndex) {
                        throw std::runtime_error("Histogram bin out of range");
                    }
                    int theEntry = thisHist + cellIndex;
                    Cache::Manager::Get()->GetHistogramsCache()
                    .SetEventIndex(eventIndex+j+1,theEntry);
                }
            }
        }
    }

    if (histCells != nextHist) {
        throw std::runtime_error("Histogram cells are missing");
    }

    return true;
}

bool Cache::Manager::Fill() {
    Cache::Manager* cache = Cache::Manager::Get();
    if (!cache) return false;
    for (auto& par : Cache::Manager::ParameterMap ) {
        if (par.first->GetParameterFunctionType()==kPolynomialCosth)
        {
            std::vector<std::vector<double>> pol_coeff = ((PolynomialCosth*)par.first->GetParameterFunction())->pol_coeff;
            int idx = par.second;
            for (int i=0;i<pol_coeff.size();i++)
            {
                for (int j=0;j<pol_coeff[i].size();j++)
                {
                    cache->GetParameterCache().SetParameter(idx++, pol_coeff[i][j]);
                }
            }
        }
        else
        {
            for (int j = 0; j < par.first->GetNpar(); j++)
                cache->GetParameterCache().SetParameter(par.second+j, par.first->GetParValue(j));
        }
    }
    cache->GetWeightsCache().Apply();
    cache->GetHistogramsCache().Apply();

#ifdef CACHE_MANAGER_SLOW_VALIDATION
#warning CACHE_MANAGER_SLOW_VALIDATION in Cache::Manager::Fill()
    // Returning false means that the event weights will also be calculated
    // using the CPU.
    return false;
#endif
    return true;
}

int Cache::Manager::ParameterIndex(const AnaFitParameters* fp) {
    std::map<const AnaFitParameters*,int>::iterator parMapIt
        = Cache::Manager::ParameterMap.find(fp);
    if (parMapIt == Cache::Manager::ParameterMap.end()) return -1;
    return parMapIt->second;
}

// An MIT Style License

// Copyright (c) 2022 Clark McGrew

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// Local Variables:
// mode:c++
// c-basic-offset:4
// compile-command:"$(git rev-parse --show-toplevel)/cmake/gundam-build.sh"
// End:
