"""
Model Deployment Package: Aneurysm Rupture Risk Prediction
Simple prediction function for clinical use
"""
import pandas as pd
import numpy as np
import pickle
import json
from pathlib import Path

class AneurysmRupturePredictor:
    """
    Clinical prediction tool for aneurysm rupture risk
    
    Usage:
        predictor = AneurysmRupturePredictor()
        risk_score, risk_category = predictor.predict(
            DE_Signature_Score=1.5,
            EC_Dysfunction_Score=0.8,
            ECM_Remodeling_Score=0.6
        )
    """
    
    def __init__(self, model_path='results/phase5_optimized_model.pkl',
                 scaler_path='results/phase5_optimized_scaler.pkl',
                 summary_path='results/phase5_best_model_summary.json'):
        """
        Initialize the predictor with trained model
        
        Args:
            model_path: Path to trained model pickle file
            scaler_path: Path to feature scaler pickle file
            summary_path: Path to model summary JSON file
        """
        # Load model and scaler
        with open(model_path, 'rb') as f:
            self.model = pickle.load(f)
        
        with open(scaler_path, 'rb') as f:
            self.scaler = pickle.load(f)
        
        with open(summary_path, 'r') as f:
            self.summary = json.load(f)
        
        self.features = self.summary['features_used']
        self.optimal_threshold = 0.605  # From phase 7 evaluation
    
    def predict(self, DE_Signature_Score=None, EC_Dysfunction_Score=None, 
                ECM_Remodeling_Score=None, **kwargs):
        """
        Predict aneurysm rupture risk
        
        Args:
            DE_Signature_Score: Differential expression signature score
            EC_Dysfunction_Score: Endothelial cell dysfunction score
            ECM_Remodeling_Score: Extracellular matrix remodeling score
            
        Returns:
            risk_score: Rupture probability (0-1)
            risk_category: 'Low', 'Medium', or 'High'
        """
        # Create feature array
        features_dict = {
            'DE_Signature_Score': DE_Signature_Score,
            'EC_Dysfunction_Score': EC_Dysfunction_Score,
            'ECM_Remodeling_Score': ECM_Remodeling_Score
        }
        
        # Fill missing values with 0
        for feature in self.features:
            if features_dict.get(feature) is None:
                features_dict[feature] = 0
        
        # Extract features in correct order
        X = np.array([[features_dict[f] for f in self.features]])
        
        # Scale features
        X_scaled = self.scaler.transform(X)
        
        # Predict probability
        risk_score = self.model.predict_proba(X_scaled)[0, 1]
        
        # Categorize risk
        if risk_score < 0.33:
            risk_category = 'Low'
        elif risk_score < 0.67:
            risk_category = 'Medium'
        else:
            risk_category = 'High'
        
        return float(risk_score), risk_category
    
    def predict_batch(self, features_df):
        """
        Predict rupture risk for multiple samples
        
        Args:
            features_df: DataFrame with columns matching feature names
            
        Returns:
            DataFrame with 'Risk_Score' and 'Risk_Category' columns
        """
        # Ensure all features are present
        for feature in self.features:
            if feature not in features_df.columns:
                features_df[feature] = 0
        
        # Extract and scale features
        X = features_df[self.features].fillna(0).values
        X_scaled = self.scaler.transform(X)
        
        # Predict probabilities
        risk_scores = self.model.predict_proba(X_scaled)[:, 1]
        
        # Categorize risks
        risk_categories = pd.cut(risk_scores, bins=[0, 0.33, 0.67, 1.0],
                                labels=['Low', 'Medium', 'High'])
        
        # Return results
        results = pd.DataFrame({
            'Risk_Score': risk_scores,
            'Risk_Category': risk_categories
        })
        
        return results
    
    def get_feature_importance(self):
        """
        Get feature importance information
        
        Returns:
            Dictionary with feature importance details
        """
        return {
            'features': self.features,
            'importance_percent': {
                'DE_Signature_Score': 56.5,
                'EC_Dysfunction_Score': 35.6,
                'ECM_Remodeling_Score': 7.8
            },
            'model_type': self.summary['best_model'],
            'optimal_threshold': self.optimal_threshold
        }


def load_nomogram(model_path='results/phase6_nomogram_formula.json'):
    """
    Load nomogram formula for combined prediction
    
    Returns:
        Nomogram information dictionary
    """
    with open(model_path, 'r') as f:
        return json.load(f)


def predict_with_nomogram(optimized_model_score, exosomal_scores=None):
    """
    Predict using combined nomogram (optimized model + exosomal biomarkers)
    
    Args:
        optimized_model_score: Score from optimized model (0-1)
        exosomal_scores: Dict with exosomal biomarker scores (optional)
        
    Returns:
        Nomogram score (0-1)
    """
    nomogram = load_nomogram()
    
    # Default exosomal scores if not provided
    if exosomal_scores is None:
        exosomal_scores = {
            'exosomal_1': 0.5,
            'exosomal_2': 0.5,
            'exosomal_3': 0.5
        }
    
    # Calculate nomogram score
    nomogram_score = (
        0.70 * optimized_model_score +
        0.15 * exosomal_scores.get('exosomal_1', 0.5) +
        0.10 * exosomal_scores.get('exosomal_2', 0.5) +
        0.05 * exosomal_scores.get('exosomal_3', 0.5)
    )
    
    return nomogram_score


# Example usage
if __name__ == '__main__':
    print("=" * 80)
    print("ANEURYSM RUPTURE RISK PREDICTOR - EXAMPLE USAGE")
    print("=" * 80)
    
    # Initialize predictor
    predictor = AneurysmRupturePredictor()
    
    # Example 1: Single prediction
    print("\nExample 1: Single Patient Prediction")
    print("-" * 80)
    
    risk_score, risk_category = predictor.predict(
        DE_Signature_Score=1.5,
        EC_Dysfunction_Score=0.8,
        ECM_Remodeling_Score=0.6
    )
    
    print(f"DE_Signature_Score: 1.5")
    print(f"EC_Dysfunction_Score: 0.8")
    print(f"ECM_Remodeling_Score: 0.6")
    print(f"\nRisk Score: {risk_score:.3f}")
    print(f"Risk Category: {risk_category}")
    
    # Example 2: Batch prediction
    print("\n\nExample 2: Batch Prediction")
    print("-" * 80)
    
    sample_features = pd.DataFrame({
        'DE_Signature_Score': [1.5, 0.2, -0.5],
        'EC_Dysfunction_Score': [0.8, 0.1, -0.3],
        'ECM_Remodeling_Score': [0.6, -0.2, -0.4]
    })
    
    results = predictor.predict_batch(sample_features)
    print(sample_features)
    print("\nPredictions:")
    print(results)
    
    # Example 3: Feature importance
    print("\n\nExample 3: Feature Importance")
    print("-" * 80)
    
    importance = predictor.get_feature_importance()
    print(f"Model Type: {importance['model_type']}")
    print(f"Optimal Threshold: {importance['optimal_threshold']}")
    print(f"\nFeature Importance:")
    for feature, imp in importance['importance_percent'].items():
        print(f"  {feature}: {imp:.1f}%")
    
    # Example 4: Nomogram prediction
    print("\n\nExample 4: Nomogram Prediction")
    print("-" * 80)
    
    optimized_score = risk_score
    nomogram_score = predict_with_nomogram(optimized_score)
    
    print(f"Optimized Model Score: {optimized_score:.3f}")
    print(f"Nomogram Score: {nomogram_score:.3f}")
    
    print("\n" + "=" * 80)
    print("PREDICTION COMPLETE")
    print("=" * 80)
