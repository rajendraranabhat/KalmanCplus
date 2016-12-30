//============================================================================
// Name        : GuogoV2.1.cpp
// Author      : Rajendra
// Version     :
// Copyright   : Your copyright notice
// Description : Guoguo in C++, Ansi-style
//============================================================================

#include <iostream>
#include "matrix.h"
#include "DATA.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm> //for max,min
#include <stdlib.h> //for abs
#include <iomanip>  //for setw
#include <Eigen/Dense>
#include<Eigen/SVD>
#include<Eigen/LU>

using namespace std;

class Guoguo{

	public:

	template<typename T>
	Matrix<T> getRanging(){
			cout<<"GuogoV2.1 "<<endl;
				Matrix<double> poskalmanx(1,2,0); //Zero matrix
				Matrix<double> poskalman(1,2,1);  //Ones matrix
				poskalman = poskalman*10;

				DATA ranging;

				//printMatrix(poskalman);
				double stations_arr[][8] = {
						{755.015, 755.015, 635.635, 516.255, 396.875, 456.565, 441.523, 632.887},
						{17.780, 167.005, 226.695, 196.850, 107.315, 17.780, 210.078,  3.810},
						{69.850, 69.850, 69.850, 69.850, 69.850, 69.850, 0,   0   }
				};
				Matrix<double> stations = ranging.getStations();


				Matrix<double> Ranging = ranging.getRangingData();
				int row = Ranging.get_rows(); //113
				bool isPrint = false;

				int Numstation = Ranging.get_cols(); //6

				Matrix<double> kalmanx(2,Numstation,0);
				Matrix<double> kalmanP(2,Numstation,1);
				kalmanP = kalmanP*10;

				Matrix<double> kalmanzeroN(1,Numstation,0);
				Matrix<double> bufferFit(1,Numstation,0);
				Matrix<double> bufferVel(1,Numstation,0);
				Matrix<double> sid(1,Numstation,0);

				double Axis_arr[] = {stations.minRow(0), stations.maxRow(0), stations.minRow(1), stations.maxRow(1)};
				Matrix<double> Axis(1,4,0);
				for(int i=0;i<4;i++){
					Axis(0,i) = Axis_arr[i];
				}

				int sid_array[] = {6,1,2,3,4,5}; //stationid
				for(int i=0;i<6;i++){
					sid(0,i) = sid_array[i];
				}

				//row = 7;

				Matrix<double> rangingV(3,Numstation,0);
				Matrix<double> rangingVFit(1,Numstation,0);
				Matrix<double> rangingVInit(1,Numstation,0);
				Matrix<double> rangingVNew(3,Numstation,0);
				Matrix<double> zeros_mat(1,Numstation,0);//Just zero matrix
				Matrix<double> target3(3,1,0);
				Matrix<double> targbuffer(row,3,0);
				Matrix<double> targetbuffer(row,2,0);

				Matrix<double> xe(4,1,0);
				Matrix<double> PP = Matrix<double>::eye(4);

			    /**********************************/
				//Inside MAIN loop
				for (int numpos=0;numpos<row;numpos++){ //row=113

					//Setting data for rangingV
					for(int j=0;j<Numstation;j++){
						rangingV(0,j) = Ranging(numpos,j);
						rangingV(1,j) = (Ranging(numpos,j)!=0);
						rangingV(2,j) = 0;
					}

					if(numpos == 0){
						rangingVFit = Matrix<double>::times(rangingV.getRow(0),rangingV.getRow(1));
						rangingVInit = rangingVFit;

						for(int i=0;i<Numstation;i++){
							rangingVNew(0,i) = rangingV(0,i);
							rangingVNew(1,i) = rangingV(1,i);
							rangingVNew(2,i) = zeros_mat(0,i); //Adding 0
						}
						//Initialize target3
						for(int i=0;i<3;i++){
							target3(0,i)=0;
						}

					}//end if
					else{
						//Second time again..
						for(int i=0;i<Numstation;i++){
							rangingV(0,i) = rangingV(0,i)*1;
							rangingV(1,i) = rangingV(1,i)*1;
							rangingV(2,i) = rangingV(0,i)-rangingVInit(rangingVInit.get_rows()-1,i);
						}

						Matrix<double> tmp_rangingV = Matrix<double>::times(rangingV.getRow(0),rangingV.getRow(1));
						rangingVInit.appendMatrix(tmp_rangingV);

						//Call ProcessTBL3 function
						rangingVNew = ProcessTBL3(numpos, Numstation, rangingV, rangingVFit.matrixRowRange(max(numpos-20,0),numpos),kalmanx,kalmanP,kalmanzeroN,bufferFit,bufferVel);

						Matrix<double> tmp_rangingVFit = Matrix<double>::times(rangingVNew.getRow(0),rangingVNew.getRow(1));
						rangingVFit.appendMatrix(tmp_rangingVFit);

					}

					/*** START LS*********************/
					int ref = 0;
					Matrix<double> targ = LS3D(rangingVNew,ref,sid,stations,target3);

					targbuffer(numpos,0) = targ(0,0);
					targbuffer(numpos,1) = targ(1,0);
					targbuffer(numpos,2) = targ(2,0);

					PositionKF3(numpos,targ,target3,xe,PP,Axis);

					targetbuffer(numpos,0) = target3(0,0);
					targetbuffer(numpos,1) = target3(1,0);

				}//end for loop


				//Eigen::MatrixXd stations_ = Matrix2Eigen<double>(stations);
				//Eigen::MatrixXd targetbuffer_ = Matrix2Eigen<double>(targetbuffer);

			return targetbuffer;
		}//end getRanging


		template<typename T>
		Matrix<T> LS3D(Matrix<T> rangingVNew, int ref, Matrix<T> sid, Matrix<T> stations, Matrix<T> target){

			int i=0,j=0,k=0;

			//Lets convert matrix to eigen type first..
			Eigen::MatrixXd rangingVNew_eig = Matrix2Eigen<double>(rangingVNew);
			Eigen::MatrixXd sid_eig = Matrix2Eigen<double>(sid);
			Eigen::MatrixXd stations_eig = Matrix2Eigen<double>(stations);
			Eigen::MatrixXd target_eig = Matrix2Eigen<double>(target);

			Eigen::MatrixXd Station(stations_eig.rows(),sid_eig.cols());
			Station.fill(0);

			//Assign corresponding station id column to Station.
			for(j=0;j<sid_eig.cols();j++){
				int sid_ = sid_eig(0,j);
				Station.col(j) = stations_eig.col(sid_-1);
			}

			Eigen::VectorXd rnew = rangingVNew_eig.row(0).cwiseProduct(rangingVNew_eig.row(1));

			//Trying to remove 0 from row. Basically trying to remove the column at all corresponding to 0 data.
			//rnew(5) = 0;
			for(int i=0;i<rnew.size();i++){
				if(rnew(i)==0){
					removeElementOfVector(rnew,i);
					removeColumnOfMatrix(Station,i);
				}
			}

			int Numstation = rnew.size();

			Eigen::MatrixXd targ_eig = target_eig;

			Matrix<double> targ = Eigen2Matrix<double>(targ_eig);

			if(Numstation < 4){

				return targ;
			}

			Eigen::MatrixXd AA(Numstation-1,4);
			Eigen::MatrixXd pp(Numstation-1,1);

			AA.fill(0);
			pp.fill(0);

			//Make 3rd row to 0 vector.
			for(j=0;j<Numstation;j++){
				Station(2,j) = 0;
			}

			double kref = pow(Station(0,ref),2) + pow(Station(1,ref),2) + pow(Station(2,ref),2);
			int id = 0;
			int ki = 0;
			for(i=0;i<Numstation;i++){
				if(i != ref){
					AA(id,0) = Station(0,i) - Station(0,ref);
					AA(id,1) = Station(1,i) - Station(1,ref);
					AA(id,2) = Station(2,i) - Station(2,ref);
					AA(id,3) = rnew(ref) - rnew(i);
					ki = pow(Station(0,i),2) + pow(Station(1,i),2) + pow(Station(2,i),2);
					pp(id,0) = ki - kref - pow( (rnew(i) - rnew(ref)), 2);
					id ++;
				}
			}

			//Use Eigen pseudoinverse to calculate pinv.
			Eigen::MatrixXd tmp_eig = AA.adjoint()*AA;

			targ_eig = pseudoInverse(tmp_eig)*AA.adjoint()*pp/2;

			targ = Eigen2Matrix<double>(targ_eig);

			return targ;
		}


		template<typename T>
		Matrix<T> ProcessTBL3(int numpos, int Numstation, Matrix<T> rangingVNew,Matrix<T> rangingVFit, Matrix<T>& kalmanx, Matrix<T>& kalmanP, Matrix<T>& kalmanzeroN,
				Matrix<T>& bufferFit, Matrix<T>& bufferVel){

			Matrix<double> rangingV = rangingVNew;

			//Converting it to Eigen matrix..
			Eigen::MatrixXd rangingVNew_eig = Matrix2Eigen<double>(rangingVNew);
			Eigen::MatrixXd rangingVFit_eig = Matrix2Eigen<double>(rangingVFit);
			Eigen::MatrixXd kalmanx_eig = Matrix2Eigen<double>(kalmanx);
			Eigen::MatrixXd kalmanP_eig = Matrix2Eigen<double>(kalmanP);
			Eigen::MatrixXd kalmanzeroN_eig = Matrix2Eigen<double>(kalmanzeroN);
			Eigen::MatrixXd bufferFit_eig = Matrix2Eigen<double>(bufferFit);
			Eigen::MatrixXd bufferVel_eig = Matrix2Eigen<double>(bufferVel);

			Eigen::MatrixXd meandelta_eig(1,Numstation);

			Eigen::MatrixXd rangingV_eig = rangingVNew_eig;

			//cout<<rangingV_eig;

			//cout<<"\n************\n";


			if(numpos < 4){
				//Do nothing
				return rangingV;
			}

			if(numpos == 4){

				int rangingVFit_eig_row = rangingVFit_eig.rows();
				rangingVFit_eig.conservativeResize(rangingVFit_eig_row+1,Eigen::NoChange);//Adding 1 extra row.
				rangingVFit_eig.row(rangingVFit_eig_row) = rangingV_eig.row(0).cwiseProduct(rangingV_eig.row(1));

				for(int j=0;j<Numstation;j++){
					double sumdelta = 0;
					int deltanum = 4;

					for(int i=1;i<5;i++){
						if(rangingVFit_eig(i,j)==0 || rangingVFit_eig(i-1,j)==0){
							deltanum = deltanum-1;
						}else{
							sumdelta = sumdelta + rangingVFit_eig(i,j)-rangingVFit_eig(i-1,j);
						}
					}
					meandelta_eig(j) = (double)sumdelta/deltanum;
				}

				kalmanx_eig.row(0) = rangingVFit_eig.row(numpos);
				kalmanx_eig.row(1) = meandelta_eig;
				bufferFit_eig = rangingVFit_eig.row(numpos);

				//Change back to regular matrix
				kalmanx = Eigen2Matrix<double>(kalmanx_eig);
				bufferFit = Eigen2Matrix<double>(bufferFit_eig);
				bufferVel = Eigen2Matrix<double>(bufferVel_eig);

				return rangingV;
			  }//End for loop

			//printMatrix(rangingVNew_eig,"rangingVNew_eig");
			Eigen::VectorXd dut1_eig(2);
			dut1_eig.fill(0); //Initialize vector with 0.
			//For numpos > 4
			for(int ii=0;ii<Numstation;ii++){
				dut1_eig(0) = rangingV_eig(0,ii)*rangingV_eig(1,ii);
				dut1_eig(1) = rangingVNew_eig(2,ii);

				rangKF3(ii, dut1_eig, kalmanx_eig, kalmanP_eig, kalmanzeroN_eig, rangingVFit_eig, bufferFit_eig, bufferVel_eig);

				rangingV_eig(0,ii) = dut1_eig(0);
				rangingV_eig(1,ii) = (dut1_eig(0) != 0);
				rangingV_eig(2,ii) = dut1_eig(1);

			}

			//Change back to regular matrix
			kalmanx = Eigen2Matrix<double>(kalmanx_eig);
			kalmanP = Eigen2Matrix<double>(kalmanP_eig);
			kalmanzeroN = Eigen2Matrix<double>(kalmanzeroN_eig);
			bufferFit = Eigen2Matrix<double>(bufferFit_eig);
			bufferVel = Eigen2Matrix<double>(bufferVel_eig);
			rangingV  = Eigen2Matrix<double>(rangingV_eig);

			return rangingV;
		}


		void rangKF3(int idx, Eigen::VectorXd& dut1_eig, Eigen::MatrixXd& kalmanx_eig, Eigen::MatrixXd& kalmanP_eig, Eigen::MatrixXd& kalmanzeroN_eig, Eigen::MatrixXd& rangingVFit_eig, Eigen::MatrixXd& bufferFit_eig, Eigen::MatrixXd& bufferVel_eig){

			// dut1_eig: ranging data;
			// kalmanx_eig: state for KF;
			// kalmanP_eig: initial is 10;
			// kalmanzeroN_eig: the number of continuous 0s.
			// bufferFitd: the last normal ranging value;
			// rp: previous ranging data;  rp(1)-> rp(2)-> rp(3)-> rp(4)-> rp(5)-> r;
			// rm0: previous measured data;
			// meandelta: the mean value of last 5 ranging data;


			Eigen::MatrixXd rp_mat = rangingVFit_eig.bottomRows(5);
			Eigen::VectorXd rp = rp_mat.col(idx);


			Eigen::MatrixXd sigma(2,2);
			sigma << 0.5, 0,
					 0  ,0.3;
			int maxStepDelta = 90;

			if (dut1_eig(0) == 0){
				kalmanzeroN_eig(0,idx) = kalmanzeroN_eig(0,idx) + 1;
				if(kalmanzeroN_eig(0,idx) > 3){
					dut1_eig(0) = 0;
					dut1_eig(1) = 0; //How to give the velocity here?
					return;
				}else{
					bufferFit_eig(0,idx) = rp(4)+ bufferVel_eig(0,idx);
					dut1_eig(0) = rp(4) + bufferVel_eig(0,idx);
					dut1_eig(1) = bufferVel_eig(0,idx);
					return;
				}
			}else{

				if(kalmanzeroN_eig(0,idx) > 3){
					if( abs(dut1_eig(0)-bufferFit_eig(0,idx)) > maxStepDelta*(kalmanzeroN_eig(0,idx) + 1 - 3) ){
						kalmanzeroN_eig(0,idx) = kalmanzeroN_eig(0,idx) + 1;
						dut1_eig(0) = 0;
						dut1_eig(1) = 0;
						return;
					}else{
						dut1_eig(1) = 0;
						kalmanx_eig.col(idx) = dut1_eig.head(2);
						sigma.setIdentity(2,2);
						kalmanzeroN_eig(0,idx) = 0;
					}
				}else{
					if( abs(dut1_eig(0)-rp(4)) > maxStepDelta){
						kalmanzeroN_eig(0,idx) = kalmanzeroN_eig(0,idx) + 1;
						if(kalmanzeroN_eig(0,idx) >3){
							dut1_eig(0) = 0;
							dut1_eig(1) = bufferVel_eig(0,idx);
							return;
						}else{
							dut1_eig(0) = rp(4) + bufferVel_eig(0,idx);
							dut1_eig(1) = bufferVel_eig(0,idx);
							return;
						}
					}else{
						if(kalmanzeroN_eig(0,idx) != 0){
							dut1_eig(1) = bufferVel_eig(0,idx);
						}
						kalmanzeroN_eig(0,idx) = 0;
					}
				}
			}//end else

			//kalman filter
			int NumSensors = 8;
			Eigen::MatrixXd F(2,2);
			Eigen::MatrixXd Q(2,2);
			Eigen::MatrixXd R(2,2);
			Eigen::MatrixXd C(2,2);
			Eigen::MatrixXd P(2,2);

			F << 1, 0.078*NumSensors,
				 0, 			   1;

			R << 50,  0,
				 20,  0;

			Q = R;

			C << 1, 0,
				 0, 1;

			P << kalmanP_eig(0,idx), 0,
				 0			       , kalmanP_eig(1,idx);


			P = F*P*F + Q;

			if( (C*P*C+R).determinant() <= 1.0e-16 ){
				kalmanx_eig.col(idx) = dut1_eig.head(2);
			}else{
				//Eigen::MatrixXd K = P*C.adjoint()*(C*P*C.adjoint()+R).inverse();
				Eigen::MatrixXd K = sigma*P*C*(C*P*C+R).inverse();

				kalmanx_eig.col(idx) = F*kalmanx_eig.col(idx) + K*(dut1_eig.head(2)-C*F*kalmanx_eig.col(idx));

				P = (Eigen::MatrixXd::Ones(2,2)-K*C)*P;

				dut1_eig.head(2) = C*kalmanx_eig.col(idx);
			}

			bufferFit_eig(0,idx) = dut1_eig(0);

			bufferVel_eig(0,idx) = dut1_eig(1);

			Eigen::MatrixXd P_tmp(2,1);
			P_tmp <<P(0,0),
					P(1,1);
			P = P_tmp;
			kalmanP_eig.col(idx) = P;
		}


		template<typename T>
		void PositionKF3(int nn, Matrix<T> targ, Matrix<T>& target, Matrix<T>& xe, Matrix<T>& P, Matrix<T> Axis){
			double TT = 0.468;
			Eigen::MatrixXd F(4,4);
			Eigen::MatrixXd Q(4,4);
			Eigen::MatrixXd R(2,2);
			Eigen::MatrixXd C(2,4);

			F << 1, 0, TT, 0,
				 0, 1, 0, TT,
				 0, 0, 1, 0,
				 0, 0, 0, 1;

			Q << pow(TT,3)/2,    pow(TT,2)/2,   0,             0,
					0,		     0,             pow(TT,3)/2,   pow(TT,2)/2,
					pow(TT,2)/2, TT,             0,            0,
					0,           0,              pow(TT,2)/2,  TT;

			R << 1,  0,
				 0,  1;

			C << 1, 0, 0, 0,
				 0, 1, 0, 0;

			double sigma = 0.5;

			//Lets convert matrix to eigen type
			Eigen::MatrixXd targ_eig = Matrix2Eigen<double>(targ);
			Eigen::MatrixXd target_eig = Matrix2Eigen<double>(target);
			Eigen::MatrixXd xe_eig = Matrix2Eigen<double>(xe);
			Eigen::MatrixXd P_eig = Matrix2Eigen<double>(P);
			Eigen::MatrixXd Axis_eig = Matrix2Eigen<double>(Axis);

			if(nn == 0){
				//target_eig.resize(targ_eig.rows(),targ_eig.cols());
				target_eig.topRows(2) = targ_eig.topRows(2);
				xe_eig.fill(0); //Fill the matrix with 0
				xe_eig.topRows(2) = target_eig.topRows(2);
				P_eig = F*P_eig*F.adjoint() + Q;

				//Conversion back..
				target = Eigen2Matrix<double>(target_eig);
				P 	   = Eigen2Matrix<double>(P_eig);
				xe	   = Eigen2Matrix<double>(xe_eig);
				return;
			}

			//cout<<"\n************\n";

			if( (targ_eig(0) <= Axis_eig(0))
					|| (targ_eig(0) >= Axis_eig(1))
					||(targ_eig(1) <= Axis_eig(2))
					||(targ_eig(1) >= Axis_eig(3))){
				return;
			}

			double rad = sqrt( pow((targ_eig(0)-target_eig(0)),2) + pow((targ_eig(1)-target_eig(1)),2));

			if(rad > 100){
				target_eig(0) = target_eig(0)*0.8 + targ_eig(0)*0.2;
				target_eig(1) = target_eig(1)*0.8 + targ_eig(1)*0.2;
				sigma = 0.1;
			}

			if(rad<5){
				target_eig(0) = target_eig(0)*0.5 + targ_eig(0)*0.5;
				target_eig(1) = target_eig(1)*0.5 + targ_eig(1)*0.5;
			}

			P_eig = F*P_eig*F.adjoint() + Q;

			double dett = (C*P_eig*C.adjoint()+R).determinant();
			if( (C*P_eig*C.adjoint()+R).determinant() <= 1.0e-30 ){
				xe_eig.topRows(2) = targ_eig.topRows(2)*0.5 + target_eig.topRows(2)*0.5;
			}else{
				Eigen::MatrixXd K_eig = P_eig*C.adjoint()*(C*P_eig*C.adjoint()+R).inverse();
				xe_eig = xe_eig + K_eig*(target_eig.topRows(2) - C*xe_eig)*sigma;
				//P_eig = (Eigen::MatrixXd::Identity(4,4)-K_eig*C)*P_eig;
				P_eig = (Eigen::MatrixXd::Ones(4,4)-K_eig*C)*P_eig;
			}

			target_eig.topRows(2) = C*xe_eig;
			target_eig(2) = targ_eig(2);

			//Conversion back..
			target = Eigen2Matrix<double>(target_eig);
			P 	   = Eigen2Matrix<double>(P_eig);
			xe	   = Eigen2Matrix<double>(xe_eig);
		}


		//For calculating pseudo inverse..
		template<typename T>
		//T pseudoInverse(const T &a, double epsilon = std::numeric_limits<double>::epsilon())
		T pseudoInverse(const T &a)
		{
			double epsilon = std::numeric_limits<double>::epsilon();
			Eigen::JacobiSVD< T > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
			//double epsilon = 0.0001;
			double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
			return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
		}

		void removeElementOfVector(Eigen::VectorXd& vector_, unsigned int idx){
			int n = vector_.size()-idx-1;
			vector_.segment(idx,n) = vector_.tail(n);
			vector_.conservativeResize(vector_.size()-1);
		}

		void removeRowOfMatrix(Eigen::MatrixXd& matrix, unsigned int rowToRemove)
		{
		    unsigned int numRows = matrix.rows()-1;
		    unsigned int numCols = matrix.cols();

		    if( rowToRemove < numRows )
		        matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = matrix.block(rowToRemove+1,0,numRows-rowToRemove,numCols);

		    matrix.conservativeResize(numRows,numCols);
		}

		void removeColumnOfMatrix(Eigen::MatrixXd& matrix, unsigned int colToRemove)
		{
		    unsigned int numRows = matrix.rows();
		    unsigned int numCols = matrix.cols()-1;

		    if( colToRemove < numCols )
		        matrix.block(0,colToRemove,numRows,numCols-colToRemove) = matrix.block(0,colToRemove+1,numRows,numCols-colToRemove);

		    matrix.conservativeResize(numRows,numCols);
		}

		template<typename T>
		Eigen::MatrixXd Matrix2Eigen(Matrix<double> tmp_targ){
			Eigen::MatrixXd eigen_mat(tmp_targ.get_rows(),tmp_targ.get_cols());
				for(int i=0;i<tmp_targ.get_rows();i++){
					for(int j=0;j<tmp_targ.get_cols();j++){
						eigen_mat(i,j) = tmp_targ(i,j);
				}
			}
			return eigen_mat;
		}

		template<typename T>
		Matrix<T> Eigen2Matrix(Eigen::MatrixXd eig_mat){
			Matrix<T> matt(eig_mat.rows(),eig_mat.cols(),0);

				for(int i=0;i<eig_mat.rows();i++){
					for(int j=0;j<eig_mat.cols();j++){
						matt(i,j) = eig_mat(i,j);
				}
			}
			return matt;
		}

		void printMatrix(Matrix<double> mat3,string varName){

			cout<<"\n*********************************"<<endl;
			cout<<varName<<endl;

			for(int i=0;i<mat3.get_rows();i++){
					for(int j=0;j<mat3.get_cols();j++){
						cout<<mat3(i,j)<<", ";
					}
					cout<<endl;
			}
			cout<<endl;
		}

		void printMatrix(Eigen::MatrixXd mat3,string varName){

			cout<<"\n*********************************"<<endl;
			cout<<varName<<endl;
			cout<<mat3;
			cout<<endl;
		}

		void printVector(Eigen::VectorXd vect3,string varName){

			cout<<"\n*********************************"<<endl;
			cout<<varName<<endl;
			cout<<vect3;
			cout<<endl;
		}

		template<typename T>
		void printArray(T *array,int rows,int cols){
			for(int i=0;i<rows;i++){
				for(int j = 0;j<cols;j++){
					cout<<setw(10)<<array[i][j]<<" ";
				}
				cout<<endl;
			}
			cout<<endl;
		}

};


//main function to test...
int main () {
	Guoguo guoguo;
	Matrix<double> ranging = guoguo.getRanging<double>();
	guoguo.printMatrix(guoguo.Matrix2Eigen<double>(ranging).adjoint(),"ranging");
	return 0;
}
























