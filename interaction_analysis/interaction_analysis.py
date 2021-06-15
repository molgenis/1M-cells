#!/usr/bin/env python3

"""
File:         test_interaction_eqtl.py
Created:      2021/01/11
Last Changed:
Author:       M.Vochteloo

Copyright (C) 2020 M.Vochteloo
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard imports.
from __future__ import print_function
from pathlib import Path
import argparse
import os

# Third party imports.
import numpy as np
import pandas as pd
import scipy.stats as stats
import statsmodels.api as sm

# Local application imports.

# Metadata
__program__ = "Test Interaction eQTL"
__author__ = "Martijn Vochteloo"
__maintainer__ = "Martijn Vochteloo"
__email__ = "m.vochteloo@rug.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)

"""
Syntax:
./test_interaction_eqtl.py -g ../../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/create_matrices/genotype_table.txt.gz -e ../../2020-10-12-deconvolution_gav/matrix_preparation/cortex_eur_cis/create_matrices/expression_table.txt.gz -c ../data/cell_fractions_combined.txt
"""


class main():
    def __init__(self):
        # Get the command line arguments.
        arguments = self.create_argument_parser()
        self.geno_path = getattr(arguments, 'genotype')
        self.expr_path = getattr(arguments, 'expression')
        self.covs_path = getattr(arguments, 'covariates')
        self.out_path = getattr(arguments, 'output')

        # Set variables.
        self.outdir = os.path.join(str(Path(__file__).parent.parent), 'test_interaction_eqtl')
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir)

    @staticmethod
    def create_argument_parser():
        parser = argparse.ArgumentParser(prog=__program__,
                                         description=__description__)

        # Add optional arguments.
        parser.add_argument("-v",
                            "--version",
                            action="version",
                            version="{} {}".format(__program__,
                                                   __version__),
                            help="show program's version number and exit")
        parser.add_argument("-g",
                            "--genotype",
                            type=str,
                            required=True,
                            help="The path to the genotype matrix")
        parser.add_argument("-e",
                            "--expression",
                            type=str,
                            required=True,
                            help="The path to the expression matrix")
        parser.add_argument("-c",
                            "--covariates",
                            type=str,
                            required=True,
                            help="The path to the covariates "
                                 "matrix")
        parser.add_argument("-o",
                            "--output",
                            type=str,
                            required=True,
                            help="The prepend to the output "
                                 "file")

        return parser.parse_args()

    def start(self):
        self.print_arguments()
        print("### STEP 1 ###")
        print("Loading input files.")
        geno_df = self.load_file(self.geno_path)
        expr_df = self.load_file(self.expr_path)
        covs_df = self.load_file(self.covs_path)

        self.validate(geno_df, expr_df, covs_df)

        print("### STEP 2 ###")
        print("Preprocessing.")
        # Replace -1 with NaN in the genotype dataframe. This way we can
        # drop missing values.
        geno_df.replace(-1, np.nan, inplace=True)

        print("### STEP 3 ###")
        print("Analyzing.")
        results = {}
        n = geno_df.shape[0]
        for row_index in range(n):
            if (row_index == 0) or (row_index % int(n / 10) == 0) or (row_index == (n-1)):
                print("\tProcessing eQTL {}/{} "
                      "[{:.0f}%]".format(row_index + 1,
                                         n,
                                         (100 / n) * (row_index + 1)))

            # Get the missing genotype indices.
            indices = np.arange(geno_df.shape[1])
            eqtl_indices = indices[~geno_df.iloc[row_index, :].isnull().values]

            # Subset the row and present samples for this eQTL.
            genotype = geno_df.iloc[row_index, eqtl_indices].copy()
            expression = expr_df.iloc[row_index, eqtl_indices].copy()

            # Create the base model.
            base_matrix = genotype.to_frame()
            base_matrix["intercept"] = 1

            # Loop over the covariates.
            for cov_index, cov_colname in enumerate(covs_df.columns):
                # Get the covariate we are processing.
                covariate = covs_df.iloc[eqtl_indices, cov_index].copy()
                cov_name = covariate.name

                # Create the null model.
                null_matrix = base_matrix.copy()
                null_matrix = null_matrix.merge(covariate, left_index=True, right_index=True)
                n_null = null_matrix.shape[0]
                df_null, rss_null, _, _, _ = self.create_model(null_matrix,
                                                               expression)

                # Create the alternative matrix and add the interaction
                # term.
                alt_matrix = null_matrix.copy()
                inter_name = "{}_X_{}".format(genotype.name, cov_name)
                alt_matrix[inter_name] = alt_matrix[genotype.name] * alt_matrix[cov_name]

                # Create the alternative model.
                n_alt = alt_matrix.shape[0]
                df_alt, rss_alt, r2_alt, coefficients_alt, std_errors_alt = self.create_model(alt_matrix, expression, cols=[inter_name])

                del null_matrix, alt_matrix

                # Make sure the n's are identical.
                if n_null != n_alt:
                    print("\t\t\tError due to unequal n_null and n_alt")
                    continue

                # Compare the null and alternative model.
                fvalue = self.calc_f_value(rss_null, rss_alt,
                                           df_null, df_alt, n_null)
                pvalue = self.get_p_value(fvalue, df_null, df_alt, n_null)

                # Update results.
                results[cov_colname] = results.get(cov_colname, []) + [["{}_{}".format(genotype.name, expression.name),
                                                                        coefficients_alt[inter_name],
                                                                        std_errors_alt[inter_name],
                                                                        fvalue,
                                                                        pvalue,
                                                                        r2_alt]]

                del fvalue, pvalue

            del indices, eqtl_indices, genotype, expression, base_matrix

        print("### STEP 4 ###")
        print("Saving results.")
        for covariate, result in results.items():
            print("Covariate: {}".format(covariate))
            df = pd.DataFrame(result, columns=["eQTL", "coef", "std", "f-value", "p-value", "rsquared"])
            df = df.set_index('eQTL')
            print(df)
            df.to_csv(os.path.join(self.out_path, "test_interaction_eqtl_{}.txt.gz".format(covariate)), compression="gzip", sep="\t", index=True, header=True)

    @staticmethod
    def load_file(path, sep="\t", header=0, index_col=0, nrows=None):
        df = pd.read_csv(path, sep=sep, header=header, index_col=index_col,
                         nrows=nrows)
        print("\tLoaded dataframe: {} "
              "with shape: {}".format(os.path.basename(path),
                                      df.shape))
        return df

    @staticmethod
    def validate(geno_df, expr_df, covs_df):
        if not covs_df.index.equals(geno_df.columns):
            print("The genotype file columns do not match the cell "
                  "type fractions file.")
            exit()

        if not covs_df.index.equals(expr_df.columns):
            print("The expressiom file columns do not match the cell "
                  "type fractions file.")
            exit()

    @staticmethod
    def create_model(X, y, cols=None):
        """
        Method for creating a multilinear model.

        :param X: DataFrame, the matrix with rows as samples and columns as
                             dimensions.
        :param y: Series, the outcome values.

        :return df: int, the degrees of freedom of this model.
        :return ssr: float, the residual sum of squares of this fit.
        :return rsquared:
        :return coefficients:
        :return std_errors:
        """
        # Perform the Ordinary least squares fit.
        ols = sm.OLS(y.values, X)
        ssr = np.nan
        rsquared = np.nan
        try:
            ols_result = ols.fit()
            ssr = ols_result.ssr
            rsquared = ols_result.rsquared
        except np.linalg.LinAlgError as e:
            print("\t\tError: {}".format(e))
            return X.shape[1], np.nan, {col: np.nan for col in cols}, {
                col: np.nan for col in cols}

        # Extract the required info from the model.
        coefficients = {}
        std_errors = {}
        if cols is not None:
            for col in cols:
                if col in X.columns:
                    coef = ols_result.params[col]
                    std_err = ols_result.bse[col]
                else:
                    coef = np.nan
                    std_err = np.nan
                coefficients[col] = coef
                std_errors[col] = std_err

        return X.shape[1], ssr, rsquared, coefficients, std_errors

    @staticmethod
    def calc_f_value(rss1, rss2, df1, df2, n):
        """
        Method for comparing the risdual sum squared of two models using
        the F statistic.

        :param rss1: float, the residual sum of squares of the null model.
        :param rss2: float, the residual sum of squares of the alternative model.
        :param df1: int, the degrees of freedom of the null model.
        :param df2: int, the degrees of freedom of the alternative model.
        :param n: int, the number of samples in the model.
        :return : float, the p-value of the comparison.
        """
        if (rss1 == np.nan) or (rss2 == np.nan):
            return np.nan
        if df1 >= df2:
            return np.nan
        if df2 >= n:
            return np.nan
        if rss2 >= rss1:
            return 0

        return ((rss1 - rss2) / (df2 - df1)) / (rss2 / (n - df2))

    @staticmethod
    def get_p_value(f_value, df1, df2, n):
        """
        Method for getting the p-value corresponding to a F-distribution.

        :param f_value: float, the f-value.
        :param df1: int, the degrees of freedom of the null model.
        :param df2: int, the degrees of freedom of the alternative model.
        :param n: int, the number of samples in the model.
        :return : float, the p-value corresponding to the f-value.
        """
        # Low
        if f_value == np.nan:
            return np.nan
        if df1 >= df2:
            return np.nan
        if df2 >= n:
            return np.nan

        return stats.f.sf(f_value, dfn=(df2 - df1), dfd=(n - df2))

    def print_arguments(self):
        print("Arguments:")
        print("  > Genotype path: {}".format(self.geno_path))
        print("  > Expression path: {}".format(self.expr_path))
        print("  > Covariates path: {}".format(self.covs_path))
        print("")


if __name__ == '__main__':
    m = main()
    m.start()
