#ifndef utilities_hxx_
#define utilities_hxx_
//#define ITK_IO_FACTORY_REGISTER_MANAGER
//#define ITK_MANUAL_INSTANTIATION
// itk
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkBinaryThresholdImageFilter.h"
#include  "itkXorImageFilter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkSampleToHistogramFilter.h"
#include "itkListSample.h"
#include "itkHistogram.h"
#include "itkDenseFrequencyContainer2.h"
#include "itkExpectationMaximizationMixtureModelEstimator.h"
#include "itkGaussianMixtureModelComponent.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"

#include <numeric>
#include <csignal>

//namespace LAWallSegmenter
  /**
   * readImage
   */
  template< typename itkImage_t >
  typename itkImage_t::Pointer readImage(const char *fileName)
  {
    typedef itk::ImageFileReader< itkImage_t > ImageReaderType;
    typename ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(fileName);

    typename itkImage_t::Pointer image;
    
    try
      {
        reader->Update();
        image = reader->GetOutput();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cerr<< "ExceptionObject caught !" << std::endl; 
        std::cerr<< err << std::endl; 
        raise(SIGABRT);
      }

    return image;
  }


  /**
   * writeImage
   */
  template< typename itkImage_t > void writeImage(typename itkImage_t::Pointer img, const char *fileName)
  {
    typedef itk::ImageFileWriter< itkImage_t > WriterType;

    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( fileName );
    writer->SetInput(img);
    writer->UseCompressionOn();

    try
      {
        writer->Update();
      }
    catch ( itk::ExceptionObject &err)
      {
        std::cout << "ExceptionObject caught !" << std::endl; 
        std::cout << err << std::endl; 
        raise(SIGABRT);   
      }
  }

template<typename InputImageType, typename OutputImageType>
typename OutputImageType::Pointer
ResampleImage(typename InputImageType::Pointer imIn, typename InputImageType::SpacingType sampleSpacing) {

    typename InputImageType::SizeType inSize;
    typename InputImageType::SpacingType inSpacing;
    typename OutputImageType::SizeType outSize;
    typename OutputImageType::SpacingType outSpacing;

    inSize = imIn->GetLargestPossibleRegion().GetSize();
    inSpacing = imIn->GetSpacing();
    outSpacing = sampleSpacing;

    for(unsigned int i = 0; i < 3; i++) {
        outSize[i] = round(inSpacing[i]*inSize[i]/outSpacing[i]);
    }

    typedef itk::ResampleImageFilter< InputImageType, OutputImageType >    ResampleFilterType;
    typedef itk::LinearInterpolateImageFunction< InputImageType, double >  InterpolatorType;

    typename ResampleFilterType::Pointer resampler = ResampleFilterType::New();
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

    resampler->SetInterpolator( interpolator );
    resampler->SetInput( imIn );
    resampler->SetSize( outSize);
    resampler->SetOutputOrigin(  imIn->GetOrigin() );
    resampler->SetOutputSpacing(sampleSpacing);
    resampler->SetOutputDirection( imIn->GetDirection() );
    resampler->SetDefaultPixelValue( 0 );
    resampler->Update();

    return resampler->GetOutput();
}

template<typename input_image_t, typename output_image_t>
typename output_image_t::Pointer
BinThreshold(typename input_image_t::Pointer input,                            \
      typename input_image_t::PixelType lowerT,                         \
      typename input_image_t::PixelType upperT, \
      typename output_image_t::PixelType insideValue = 1,               \
      typename output_image_t::PixelType outsideValue = 0)
{
  /**
   * O(x) :=    I(x) \in [lowerT, upperT] ? insideValue : outsideValue
   */

  //tst
  //   std::cout<<lowerT<<std::endl;
  //   std::cout<<upperT<<std::endl;

  typedef itk::BinaryThresholdImageFilter<input_image_t, output_image_t> binaryThresholdImageFilter_t;

  typename binaryThresholdImageFilter_t::Pointer thlder = binaryThresholdImageFilter_t::New();
  thlder->SetInput(input);
  thlder->SetInsideValue(insideValue);
  thlder->SetOutsideValue(outsideValue);
  thlder->SetUpperThreshold(upperT);
  thlder->SetLowerThreshold(lowerT);
  thlder->Update();

  return thlder->GetOutput();
}

/**
    Compute XOR image for extracting myocardium wall
    NOTE: image resampling is needed in case the input images have different SIZEs!
*/
template<typename inputType1, typename inputType2, typename outputType>
typename outputType::Pointer
ExtractXORImage(typename inputType1::Pointer input1, typename inputType2::Pointer input2) {

        // Binarize masks
        input1 = BinThreshold<inputType1, inputType1>(input1, 1, 10000);
        input2 = BinThreshold<inputType2, inputType2>(input2,  1, 10000);

        typedef itk::XorImageFilter <outputType> XorImageFilterType;
        typename XorImageFilterType::Pointer xorFilter = XorImageFilterType::New();
        xorFilter->SetInput1(input1);
        xorFilter->SetInput2(input2);
        xorFilter->Update();

        return  xorFilter->GetOutput();
}

template<typename MaskImType>
typename MaskImType::Pointer
BinaryErode(typename MaskImType::Pointer maskImg, int radius) {

    typedef itk::BinaryBallStructuringElement<typename MaskImType::PixelType, 3> StrEleType;
    typedef itk::BinaryErodeImageFilter<MaskImType, MaskImType, StrEleType> ErodeFilterType;

    typename ErodeFilterType::Pointer  binaryErode  = ErodeFilterType::New();

    StrEleType  strEle;
    strEle.SetRadius(radius);   // structuring element
    strEle.CreateStructuringElement();

    binaryErode->SetKernel( strEle);
    binaryErode->SetInput(maskImg);
    binaryErode->SetForegroundValue(1);
    binaryErode->SetBackgroundValue(0);
    binaryErode->Update();

    return binaryErode->GetOutput();
}

template<typename MaskImType>
typename MaskImType::Pointer
BinaryDilate(typename MaskImType::Pointer maskImg, int radius) {

    typedef itk::BinaryBallStructuringElement<typename MaskImType::PixelType, 3> StrEleType;
    typedef itk::BinaryDilateImageFilter<MaskImType, MaskImType, StrEleType> DilateFilterType;

    typename DilateFilterType::Pointer binaryDilate = DilateFilterType::New();

    StrEleType  strEle;
    strEle.SetRadius(radius);   // structuring element
    strEle.CreateStructuringElement();

    binaryDilate->SetKernel( strEle);
    binaryDilate->SetInput(maskImg);
    binaryDilate->SetForegroundValue(1);
    binaryDilate->SetBackgroundValue(0);
    binaryDilate->Update();

    return binaryDilate->GetOutput();
}

template<typename MaskImType>
typename MaskImType::Pointer
BinaryClosing(typename MaskImType::Pointer maskImg, int radius) {

    maskImg = BinaryDilate<MaskImType>(maskImg, radius);
    maskImg = BinaryErode<MaskImType>(maskImg, radius);

    return maskImg;
}

template<typename MaskImType>
typename MaskImType::Pointer
BinaryOpening(typename MaskImType::Pointer maskImg, int radius) {

    maskImg = BinaryErode<MaskImType>(maskImg, radius);
    maskImg = BinaryDilate<MaskImType>(maskImg, radius);

    return maskImg;
}

/**
    1D Expectation Maximization
*/
template<typename TData, typename TPara>
void ExpectationMaximization1D(const std::vector<TData>& vec_data, std::vector<TPara>& vec_para, const unsigned int NCLASS) {

    // rough initialization
    std::vector<TPara> means(NCLASS,0.0), variance(NCLASS, 0.0), nums(NCLASS,0.0);
    for(unsigned int i = 0; i < NCLASS; i++) {
        for(unsigned int j = (unsigned int)((float)(vec_data.size())/NCLASS*i); \
            j < std::min((unsigned int)((float)(vec_data.size())/NCLASS*(i+1)), (unsigned int)(vec_data.size())); j++) {
            means[i] += vec_data[j];
            variance[i] += vec_data[j]*vec_data[j];
            nums[i]++;
        }
    }

    for(unsigned int i = 0; i < NCLASS; i++) {
        variance[i] = (nums[i]*variance[i] - means[i]*means[i])/(nums[i]*(nums[i] - 1));
        means[i] /= nums[i];
    }

    // iterations
    const unsigned int COUNT = 200;
    const float PI = 3.1415;
    std::vector<float> lambda(NCLASS, 1.0/NCLASS);

    std::vector<std::vector<TPara> > compValAll(vec_data.size(),std::vector<TPara>(NCLASS));
    TPara pro2, compSum;

    for(unsigned int k = 0; k < COUNT; k++) {

        // E-Step
        for(unsigned int j = 0; j < NCLASS; j++) {
            pro2 = 1.0/(std::sqrt(2*PI*variance[j]));
            for(unsigned int i = 0; i < vec_data.size(); i++) {
                compValAll[i][j] = lambda[j]*pro2*std::exp(-0.5*(vec_data[i] - means[j])*(vec_data[i] - means[j])/variance[j]);
            }
        }

        for(unsigned int i = 0; i < vec_data.size(); i++) {
            compSum = 0.0;
            for(unsigned int j = 0; j < NCLASS; j++) {
                compSum += compValAll[i][j];
            }
            for(unsigned int j = 0; j < NCLASS; j++) {
                compValAll[i][j] /= compSum;
            }
        }

        // M-Step
        for(unsigned int j = 0; j < NCLASS; j++) {
            lambda[j] = 0.0;
            for(unsigned int i = 0; i < vec_data.size(); i++) {
                lambda[j] += compValAll[i][j];
            }
            lambda[j] /= vec_data.size();

            means[j] = 0;
            for(unsigned int i = 0; i < vec_data.size(); i++) {
                means[j] += compValAll[i][j]*vec_data[i];
            }
            means[j] = means[j]/(vec_data.size()*lambda[j]);

            variance[j] = 0;
            for(unsigned int i = 0; i < vec_data.size(); i++) {
                variance[j] += compValAll[i][j]*(vec_data[i] - means[j])*(vec_data[i] - means[j]);
            }
            variance[j] = variance[j]/(vec_data.size()*lambda[j]);
        }
    }

    vec_para.resize(NCLASS*2);
    for(unsigned int i = 0; i < NCLASS; i++) {
        vec_para[i] = means[i];
        vec_para[NCLASS + i] = variance[i];
    }

}

template<typename MaskImType>
typename MaskImType::Pointer
CleanBinaryImage(typename MaskImType::Pointer mask, const int SIZEMIN) {

    typedef unsigned short CCPixel;
    typedef itk::Image<CCPixel, 3> CCImType;
    typedef itk::ConnectedComponentImageFilter<MaskImType, CCImType> CCFilter;
    typedef itk::RelabelComponentImageFilter<CCImType, CCImType> ReLabelFilter;

    typename CCFilter::Pointer ccFilter = CCFilter::New();
    typename ReLabelFilter::Pointer relablFilter = ReLabelFilter::New();
    //		 ccFilter->SetFullyConnected(true);
    ccFilter->SetInput(mask);		// Find Connected Component

    relablFilter->SetInput(ccFilter->GetOutput());      // Relabel Connected Component

        std::vector<float> objSize = relablFilter->GetSizeOfObjectsInPhysicalUnits();

        std::cout << "num of objects: " << objSize.size() << std::endl;
         for(unsigned int i = 0; i < objSize.size(); i++) {
             std::cout << "class " << i << ":" << objSize[i] << std::endl;
         }

     // TODO: remove small objects

     return mask;
}

/**
   Scar identification based on EM
*/
template<typename imSrcType, typename imMaskType, typename imOutputType>
typename imOutputType::Pointer
ScarSegmentation(typename imSrcType::Pointer imSrc, typename imMaskType::Pointer imWall, typename imMaskType::Pointer imLA) {

    typedef itk::Vector<typename imSrcType::PixelType, 1> MeasurementVectorType;
    typedef itk::Statistics::ListSample< MeasurementVectorType> SampleType;

    typedef itk::Statistics::Histogram< float,itk::Statistics::DenseFrequencyContainer2 > HistogramType;

    typename SampleType::Pointer wallSample = SampleType::New();
    typename SampleType::Pointer laSample = SampleType::New();

    // compute LA and Wall histograms
    typedef itk::ImageRegionConstIterator<imSrcType> ImgRegionConstInteratorSrc;
    ImgRegionConstInteratorSrc itSrc(imSrc, imSrc->GetLargestPossibleRegion());

    typedef itk::ImageRegionConstIterator<imMaskType> ImgRegionConstInteratorWall;
    ImgRegionConstInteratorWall itWall(imWall, imWall->GetLargestPossibleRegion());

    typedef itk::ImageRegionConstIterator<imMaskType> ImgRegionConstInteratorLA;
    ImgRegionConstInteratorLA itLA(imLA, imLA->GetLargestPossibleRegion());

    for(itSrc.GoToBegin(),itWall.GoToBegin(), itLA.GoToBegin();
        !itSrc.IsAtEnd(); ++itWall, ++itSrc, ++itLA) {

          if(itWall.Get() > 0) {
              wallSample->PushBack(itSrc.Get());
          }
          else if(itLA.Get() > 0) {
              laSample->PushBack(itSrc.Get());
          }
    }

    typedef itk::Statistics::SampleToHistogramFilter<SampleType, HistogramType> SampleToHistogramFilterType;
    typename SampleToHistogramFilterType::Pointer wallHist = SampleToHistogramFilterType::New();
    wallHist->SetInput(wallSample);
    typename SampleToHistogramFilterType::HistogramSizeType histogramSize(1);
    histogramSize.Fill(100);
    wallHist->SetHistogramSize(histogramSize);
    wallHist->Update();

    typename SampleToHistogramFilterType::Pointer laHist = SampleToHistogramFilterType::New();
    laHist->SetInput(laSample);
    histogramSize.Fill(100);
    laHist->SetHistogramSize(histogramSize);
    laHist->Update();

    const HistogramType* histLA = laHist->GetOutput();

    // resample histogram inside the wall
    const float PROLAMAX = 0.9;
    float binMean;
    unsigned int freqMax = 0, freq;


    for(unsigned int i = 0; i < histLA->GetSize()[0]; i++)  {
        if(histLA->GetFrequency(i) > freqMax) {
            freqMax = histLA->GetFrequency(i);
        }
    }

    const HistogramType* histWall = wallHist->GetOutput();
    std::vector<typename imSrcType::PixelType> wallResample;
    std::vector<float> wallParas(4);

    for(unsigned int i = 0; i < histWall->GetSize()[0]; i++)  {
        freq = histWall->GetFrequency(i)*(1.0 - (float)(histLA->GetFrequency(i))/freqMax*PROLAMAX);
        binMean = (histWall->GetBinMin(0, i) + histWall->GetBinMax(0, i))/2.0;
        for(unsigned int j = 0; j < freq; j++) {
            wallResample.push_back(binMean);
        }
    }

    const unsigned int NCLASS = 2;
    ExpectationMaximization1D<typename imSrcType::PixelType, float>(wallResample, wallParas, NCLASS);

    // pick the scar
    float scarTH;
    scarTH = std::max(wallParas[0], wallParas[1]);
//    std::cout << "scarMean: " << scarTH << std::endl;

    typename imOutputType::Pointer imScar = imOutputType::New();
    imScar->SetRegions(imWall->GetLargestPossibleRegion() );
    imScar->CopyInformation(imWall);
    imScar->Allocate();
    imScar->FillBuffer(0);

    typedef itk::ImageRegionIterator<imMaskType> ImgRegionInteratorScar;
    ImgRegionInteratorScar itScar(imScar, imScar->GetLargestPossibleRegion());

    for(itSrc.GoToBegin(),itWall.GoToBegin(), itScar.GoToBegin();
        !itSrc.IsAtEnd(); ++itWall, ++itSrc, ++itScar) {

        if(itWall.Get() > 0 && itSrc.Get() > scarTH) {
            itScar.Set(1);
        }
    }

    return imScar;
}



#endif
