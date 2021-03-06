Cs273aproject::Application.routes.draw do
  # The priority is based upon order of creation: first created -> highest priority.
  # See how all your routes lay out with "rake routes".

  # You can have the root of your site routed with "root"
  # root 'welcome#index'

  # Example of regular route:
  #   get 'products/:id' => 'catalog#view'

  # Example of named route that can be invoked with purchase_url(id: product.id)
  #   get 'products/:id/purchase' => 'catalog#purchase', as: :purchase

  # Example resource route (maps HTTP verbs to controller actions automatically):
  #   resources :products

  # Example resource route with options:
  #   resources :products do
  #     member do
  #       get 'short'
  #       post 'toggle'
  #     end
  #
  #     collection do
  #       get 'sold'
  #     end
  #   end

  # Example resource route with sub-resources:
  #   resources :products do
  #     resources :comments, :sales
  #     resource :seller
  #   end

  # Example resource route with more complex sub-resources:
  #   resources :products do
  #     resources :comments
  #     resources :sales do
  #       get 'recent', on: :collection
  #     end
  #   end

  # Example resource route with concerns:
  #   concern :toggleable do
  #     post 'toggle'
  #   end
  #   resources :posts, concerns: :toggleable
  #   resources :photos, concerns: :toggleable

  # Example resource route within a namespace:
  #   namespace :admin do
  #     # Directs /admin/products/* to Admin::ProductsController
  #     # (app/controllers/admin/products_controller.rb)
  #     resources :products
  #   end
  
  resources :home

  resources :admin do
    post 'run_downloader', on: :collection
  end

  resources :analysis do
    post 'stageoutput', on: :collection
    get 'stage', on: :collection
    post 'stage', on: :collection
    get 'survival', on: :collection
    post 'survival', on: :collection
    get 'get_features', on: :collection
    get 'batch_ids', on: :collection
    get 'get_prediction_targets', on: :collection
    get 'get_clinical_variables', on: :collection
    get 'get_unique_prediction_target_values', on: :collection
  end

  resources :clinical do
	  get :chart_data, on: :collection
    get 'get_clinical_variables', on: :collection
    get :get_data_sources, on: :collection
    get :table_data, on: :collection
    get :get_chart_type, on: :collection
    get :table_data_continuous, on: :collection
  end

  resources :omics do
    post :index, on: :collection
    get :get_gene_names, on: :collection
  end
  
  resources :results do
    get :progress, on: :collection
  end
  
  resources :acknowledgements
end
