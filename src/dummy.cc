/*
 * dummy.cc
 *
 * License: Artistic License, see file LICENSE.TXT or 
 *          https://opensource.org/licenses/artistic-license-1.0
 */

#include "intronmodel.hh"
#include "exonmodel.hh"
#include "igenicmodel.hh"
#include "utrmodel.hh"

void IntronModel::buildModel( const AnnoSequence* annoseq, int parIndex){}
void IntronModel::printProbabilities( int zusNumber, BaseCount *bc, const char* suffix ){}
void IntronModel::registerPars( Parameters* parameters ){}
void IntronModel::storeGCPars(int idx){}
void ExonModel::buildModel( const AnnoSequence* annoseq, int parIndex){}
void ExonModel::printProbabilities( int zusNumber, BaseCount *bc, const char* suffix ){}
void ExonModel::registerPars( Parameters* parameters ){}
void ExonModel::storeGCPars(int idx){}
void IGenicModel::buildModel( const AnnoSequence* annoseq, int parIndex){}
void IGenicModel::printProbabilities( int zusNumber, BaseCount *bc, const char* suffix ){}
void IGenicModel::registerPars( Parameters* parameters ){}
void IGenicModel::storeGCPars(int idx){}
void UtrModel::buildModel( const AnnoSequence* annoseq, int parIndex){}
void UtrModel::printProbabilities( int zusNumber, BaseCount *bc, const char* suffix ){}
void UtrModel::registerPars( Parameters* parameters ){}
void UtrModel::storeGCPars(int idx){}
